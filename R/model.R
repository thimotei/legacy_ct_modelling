params_avail_to_adjust <- function(params = "all") {
    choices <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")
    params <- match.arg(params, c(choices, "all"), several.ok = TRUE)
    if (any(params %in% "all")) {
      params <- choices
    }
    params_list <- as.list(choices)
    names(params_list) <- choices
    params_list <- purrr::map(params_list, ~ as.numeric(any(params %in% .)))
    return(params_list)
  }

  test_design <- function(formula = ~ 1, data, preds_sd = 1) {
    design <- model.matrix(formula, data = data)

    out <- list(
      design = design, preds_sd = preds_sd
    )
    return(out)
  }

subject_design <- function(formula = ~ 1, data, preds_sd = 0.1,
                             params = "all") {
  params <- params_avail_to_adjust(params)

  subjects <- extract_subjects(data)
  design <- model.matrix(formula, data = subjects)

  out <- list(
    design = design, subjects = subjects, params = params,
    preds_sd = preds_sd
  )
  return(out)
}

get_inc_period <- function(inc_mean = c(1.621, 0.0640),
                           inc_sd = c(0.418, 0.0691)) {
  list(
    inc_mean_p = inc_mean,
    inc_sd_p = inc_sd
  )
}

epict_to_stan <- function(obs,
                         ct_model = subject_design(~ 1, obs),
                         adjustment_model = test_design(~ 1, obs),
                         individual_variation = 0.2,
                         individual_correlation = 1,
                         censoring_threshold = 40, 
                         positivity_threshold = 37,
                         switch = TRUE,
                         onsets = TRUE, incubation_period = get_inc_period(),
                         likelihood = TRUE, output_loglik = FALSE) {
  obs <- data.table::copy(obs)
  obs <- obs[order(id)]
  obs[, obs := 1:.N]

  tests_per_id <- obs[, .(n = .N), by = "id"]$n

  stan_data <- list(N = obs[, .N],
                    P = length(unique(obs$id)),
                    id = obs[, id],
                    tests_per_id = tests_per_id,
                    cum_tests_per_id = cumsum(tests_per_id),
                    day_rel = obs[, t],
                    ct_value = obs$ct_value,
                    ncensored = length(obs[uncensored == 0, obs]),
                    censored = obs[uncensored == 0, obs],
                    nuncensored = length(obs[uncensored == 1, obs]),
                    uncensored = obs[uncensored == 1, obs],
                    uncensored_by_test = obs[, uncensored],
                    t_e = 0,
                    c_0 = censoring_threshold,
                    c_lod = censoring_threshold,
                    K = ifelse(switch, 5, 3),
                    ind_var_sd = individual_variation,
                    ind_var_m = as.numeric(individual_variation != 0),
                    ind_corr = as.numeric(!is.na(individual_correlation) &&
                      individual_variation != 0),
                    lkj_prior = ifelse(
                      is.na(individual_correlation), 0, individual_correlation
                    ),
                    lmean = incubation_period$inc_mean_p,
                    lsd = incubation_period$inc_sd_p,
                    likelihood = as.numeric(likelihood),
                    output_loglik = as.numeric(output_loglik),
                    preds = ncol(ct_model$design) - 1,
                    preds_sd = ct_model$preds_sd,
                    design = ct_model$design,
                    switch = as.numeric(switch),
                    adj_t_p = ct_model$params[["t_p"]],
                    adj_t_s = min(ct_model$params[["t_s"]], as.numeric(switch)),
                    adj_t_lod = ct_model$params[["t_lod"]],
                    adj_c_p = ct_model$params[["c_p"]],
                    adj_c_s = min(ct_model$params[["c_s"]], as.numeric(switch)),
                    adj_inc_mean = ct_model$params[["inc_mean"]],
                    adj_inc_sd = ct_model$params[["inc_sd"]],
                    adj_ct = as.numeric(
                      (ncol(adjustment_model$design) - 1) > 0
                    ),
                    ct_preds = ncol(adjustment_model$design) - 1,
                    ct_preds_sd = adjustment_model$preds_sd,
                    ct_design = adjustment_model$design
  )
 if (is.null(obs$onset_time) | !onsets) {
  stan_data <- c(stan_data, list(
          any_onsets = 0,
          onset_avail = rep(0, stan_data$P),
          nonsets = 1,
          ids_with_onsets = as.array(0),
          onset_time = rep(0, stan_data$P),
          onset_window = rep(0, stan_data$P)
        ))
 }else{
  onset_dt <- suppressWarnings(
    obs[,
    .(onset_time = min(onset_time, na.rm = TRUE), id), by = "id"
    ][
      is.infinite(onset_time), onset_time := NA
    ]
  )
  stan_data <- c(stan_data, list(
          any_onsets = 1,
          onset_avail = as.numeric(!is.na(onset_dt$onset_time)),
          nonsets = sum(as.numeric(!is.na(onset_dt$onset_time))),
          ids_with_onsets = onset_dt[!is.na(onset_time), id],
          onset_time = onset_dt$onset_time %>%
            tidyr::replace_na(0),
          onset_window = rep(1, stan_data$P)
        ))
 }

  return(stan_data)
}

epict_inits <- function(dt) {
  function() {
    inits <- list(
      t_inf = purrr::map_dbl(
        1:dt$P,
        ~ truncnorm::rtruncnorm(
          1, a = max(-dt$onset_time[.], 0),
          mean = max(-dt$onset_time[.] + 5, 5), sd = 1
        )
      ),
      c_0 = truncnorm::rtruncnorm(
        1, a = dt$c_lod, mean = dt$c_lod + 10, sd = 1
      ),
      c_p_mean = rnorm(1, 0, 1),
      t_p_mean = rnorm(1, 1.61, 0.5),
      t_lod_mean = rnorm(1, 2.3, 0.5),
      ind_var = abs(rnorm(dt$K, 0, dt$ind_var_sd * 0.1)),
      ind_eta  = matrix(rnorm(dt$P * dt$K, 0, 0.1), nrow = dt$K, ncol = dt$P),
      sigma = truncnorm::rtruncnorm(1, a = 0, mean = 5, sd = 0.5)
    )

    if (dt$switch > 0) {
      inits$c_s_mean <- array(rnorm(1, 0, 1))
      inits$t_s_mean <- array(rnorm(1, 1.61, 0.5))
    }

    if (dt$preds > 0) {
      if (dt$adj_t_p > 0) {
        inits$beta_t_p <- rnorm(dt$preds, 0, 0.01)
      }
      if (dt$adj_t_s > 0) {
        inits$beta_t_s <- rnorm(dt$preds, 0, 0.01)
      }
      if (dt$adj_t_lod > 0) {
        inits$beta_t_lod <- rnorm(dt$preds, 0, 0.01)
      }
      if (dt$adj_c_p > 0) {
        inits$beta_c_p <- rnorm(dt$preds,  0, 0.01)
      }
      if (dt$adj_c_s > 0) {
        inits$beta_c_s <- rnorm(dt$preds,  0, 0.01)
      }
      if (dt$adj_inc_mean > 0) {
        inits$beta_inc_mean <- rnorm(dt$preds, 0, 0.01)
      }
      if (dt$adj_inc_sd > 0) {
        inits$beta_inc_sd <- rnorm(dt$preds, 0, 0.01)
      }
    }

    if (dt$ct_preds > 0) {
      if (dt$adj_ct > 0) {
        inits$beta_ct_shift <- rnorm(dt$ct_preds, 0, 0.001)
        inits$beta_ct_scale <- rnorm(dt$ct_preds, 0, 0.001)
      }
    }

    if (dt$any_onsets == 1) {
      inits$inc_mean <- rnorm(1, dt$lmean[1], dt$lmean[2] * 0.1)
      inits$inc_sd <- truncnorm::rtruncnorm(
        1, a = 0, mean = dt$lsd[1], sd = dt$lsd[2] * 0.1
      )
    }
    return(inits)
  }
}

load_epict_model <- function() {
  mod <- cmdstan_model(
    "stan/ct_trajectory_model.stan",
    include_paths = "stan",
    stanc_options = list("O1")
 )
 return(mod)
}