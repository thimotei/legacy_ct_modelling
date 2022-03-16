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

data_to_stan <- function(input_data,
                         ct_model = subject_design(~ 1, input_data),
                         adjustment_model = test_design(~ 1, input_data),
                         likelihood = TRUE, clod = 40,
                         onsets = TRUE, correlation = 1) {
  input_data <- data.table::copy(input_data)
  input_data <- input_data[order(id)]

  tests_per_id <- input_data[, .(n = .N), by = "id"]$n

  stan_data <- list(N = input_data[, .N],
                    P = length(unique(input_data$id)),
                    id = input_data[, id],
                    tests_per_id = tests_per_id,
                    cum_tests_per_id = cumsum(tests_per_id),
                    day_rel = input_data[, t],
                    ct_value = ifelse(
                      is.na(input_data$ct_value), -99, input_data$ct_value
                    ),
                    pcr_res = input_data[, pcr_res],
                    t_e = 0,
                    c_0 = clod,
                    c_lod = clod,
                    K = 5,
                    lkj_prior = correlation,
                    lmean = get_inc_period()$inc_mean_p,
                    lsd = get_inc_period()$inc_sd_p,
                    likelihood = as.numeric(likelihood),
                    preds = ncol(ct_model$design) - 1,
                    preds_sd = ct_model$preds_sd,
                    design = ct_model$design,
                    adj_t_p = ct_model$params[["t_p"]],
                    adj_t_s = ct_model$params[["t_s"]],
                    adj_t_lod = ct_model$params[["t_lod"]],
                    adj_c_p = ct_model$params[["c_p"]],
                    adj_c_s = ct_model$params[["c_s"]],
                    adj_inc_mean = ct_model$params[["inc_mean"]],
                    adj_inc_sd = ct_model$params[["inc_sd"]],
                    adj_ct = (ncol(adjustment_model$design) - 1) > 0,
                    ct_preds = ncol(adjustment_model$design) - 1,
                    ct_preds_sd = adjustment_model$preds_sd,
                    ct_design = adjustment_model$design
  )
 if (is.null(input_data$onset_time) | !onsets) {
  stan_data <- c(stan_data, list(
          any_onsets = 0,
          onset_avail = rep(0, stan_data$P),
          onset_time = rep(0, stan_data$P)
        ))
 }else{
  onset_dt <- suppressWarnings(
    input_data[,
    .(onset_time = min(onset_time, na.rm = TRUE)), by = "id"
    ][
      is.infinite(onset_time), onset_time := NA
    ]
  )
  stan_data <- c(stan_data, list(
          any_onsets = 1,
          onset_avail = as.numeric(!is.na(onset_dt$onset_time)),
          onset_time = onset_dt$onset_time %>%
            tidyr::replace_na(0)
        ))
 }

  return(stan_data)
}

stan_inits <- function(dt) {
  function() {
    inits <- list(
      T_e = purrr::map_dbl(
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
      c_s_mean = rnorm(1, 0, 1),
      t_p_mean = rnorm(1, 1.61, 0.5),
      t_s_mean = rnorm(1, 1.61, 0.5),
      t_lod_mean = rnorm(1, 2.3, 0.5),
      ind_var = abs(rnorm(dt$K, 0, 0.1)),
      ind_eta  = matrix(rnorm(dt$P * dt$K, 0, 1), nrow = dt$K, ncol = dt$P),
      sigma = truncnorm::rtruncnorm(1, a = 0, mean = 5, sd = 0.5)
    )

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
        inits$beta_ct_shift <- rnorm(dt$ct_preds, 0, 0.01)
        inits$beta_ct_scale <- rnorm(dt$ct_preds, 0, 0.01)
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
