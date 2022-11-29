summarise_pop_pp <- function(fit) {
  draws <- fit$summary(
    variables = c(
      "t_p_mean", "t_s_mean[1]", "t_lod_mean", "c_p_mean",
      "c_s_mean[1]", "inc_mean", "inc_sd"
    )
  )
  draws <- data.table::as.data.table(draws)
  return(draws[])
}

summarise_coeff_pp <- function(fit, params, exponentiate = FALSE) {

  params <- params_avail_to_adjust(params)
  params <- names(params[purrr::map_lgl(params, ~ . == 1)])

  draws <- fit$summary(
    variables = paste0("beta_", params)
  )
  draws <- data.table::as.data.table(draws)

  if (exponentiate) {
    beta_cols <-  c("mean", "median", "sd", "mad", "q5", "q95")
    draws[, (beta_cols) := lapply(.SD, exp), .SDcols = beta_cols]
  }
  return(draws[])
}

summarise_effects <- function(draws, design, variables, exponentiate = TRUE) {
  eff_draws <- extract_coeffs(
      draws, exponentiate = exponentiate, design = design, variables
    )

    by <- "variable"
    if (!missing(design)) {
      by <- c(by, "predictor")
    }
    eff_summary <- summarise_draws(eff_draws, by = by)
    return(eff_summary)
}

summarise_adjustment <- function(draws, design) {
  eff_draws <- extract_coeffs(
    draws, exponentiate = FALSE, design = design,
    variables = c("ct_shift", "ct_scale")
  )

  eff_draws[variable %in% "ct_scale", value := exp(value)]

  by <- "variable"
  if (!missing(design)) {
    by <- c(by, "predictor")
  }
  eff_summary <- summarise_draws(eff_draws, by = by)
  return(eff_summary)
}

summarise_pp <- function(fit, obs) {
  ct_pp <- extract_posterior_predictions(fit, obs)
  ct_pp <- summarise_draws(
    ct_pp[, value := sim_ct], by = c("id", "t", "obs")
  )
  return(ct_pp)
}

summarise_draws <- function(draws, by = c("id", "time")) {

  out <- draws[,
      .(median = quantile(value, c(0.5), na.rm = TRUE),
        lo95 = quantile(value, c(0.05), na.rm = TRUE),
        lo90 = quantile(value, c(0.05), na.rm = TRUE),
        lo60 = quantile(value, c(0.20), na.rm = TRUE),
        lo30 = quantile(value, c(0.35), na.rm = TRUE),
        hi30 = quantile(value, c(0.65), na.rm = TRUE),
        hi60 = quantile(value, c(0.80), na.rm = TRUE),
        hi90 = quantile(value, c(0.95), na.rm = TRUE),
        hi95 = quantile(value, c(0.95), na.rm = TRUE)
      ),
      by = by
    ]
  return(out[])
}

summarise_ct_traj <- function(dt,
                              sum_cols = c("value",
                                           "t"),
                              lower_ct_limit,
                              censor_flag = TRUE,
                              pop_flag = TRUE) {
  
  dt_proc <- data.table::copy(dt)
  
  if(pop_flag == TRUE) { 
    sum_cols = c(sum_cols, by = "predictor")
    } 
  else {
    sum_cols = c(sum_cols, by = "id") 
  }
  
  # summarising ct trajectories
  out <- summarise_draws(dt[,
  value := ct_value
  ][, 
    ..sum_cols
  ], by = setdiff(sum_cols, "value")
  )
  
  # censor
  out <- out[
    median > 40, median := 40][
    lo90 > 40, lo90 := 40][
    hi90 > 40, hi90 := 40]
  
  # adding regressor categories
  if(pop_flag == TRUE) { 
    out <- add_regressor_categories(out)
  }
  
  return(out)
}

summarise_effect_sizes_natural <- function(draws) {
  
  draws <- extract_draws(fit)
  
  adj_params <- c("c_p", "t_p", "t_lod")
  baseline_predictors <- c("Omicron",
                           "Symptomatic",
                           "4 exposures",
                           "Age: 35-49")
  
  adj_draws <- adjust_params(draws, design = ct_model$design, onsets = FALSE) 
  
  adj_draws <- adj_draws %>% 
    update_predictor_labels()
  
  baseline_description <- "Baseline:
  Omicron,  
  4 exposures,
  symptomatic,
  Age:35—49"
  
  # adding baseline case with long description, for legend and vertical lines
  adj_draws <- add_baseline_to_draws(adj_draws,
                                     baseline_description,
                                     onsets_flag = FALSE)
  
  # adding the baseline results again for each subcategory, to plot
  # the results separately for each regressor category
  adj_draws <- add_baseline_to_draws(adj_draws, "Omicron (BA.1)", onsets_flag = F)
  adj_draws <- add_baseline_to_draws(adj_draws, "4 exposures", onsets_flag = F)
  adj_draws <- add_baseline_to_draws(adj_draws, "Symptomatic", onsets_flag = F)
  adj_draws <- add_baseline_to_draws(adj_draws, "Age: 35-49", onsets_flag = F)
  
  pop_draws <- extract_ct_params(adj_draws, by = "predictor", mean = FALSE) %>% 
    transform_to_model()
  
  pop_draws_long <- melt(pop_draws[, !"c_0"],
                         measure.vars = adj_params)
  
  effect_size_summary_natural <- pop_draws_long %>% 
    update_variable_labels() %>% 
    summarise_draws(by = c("variable", "predictor")) %>% 
    add_regressor_categories()
  
  return(effect_size_summary_natural)
}

summarise_shedding_trajectories <- function(dt) {
  
  dt_proc <- data.table::copy(dt)
  
  out <- dt_proc[, .(me = quantile(auc_prop, 0.5),
                     lo = quantile(auc_prop, 0.05),
                     hi = quantile(auc_prop, 0.95)),
                 by = c("t", "predictor", "regressor_category")]
  
  return(out)
}

summarise_cumulative_shedding_effects <- function(dt) {

  dt_proc <- copy(dt)
  
  out <- dt_proc[, .SD[which.min(abs(t - median))],
                 by = c("predictor", "regressor_category")]
  
  return(out)
}

summarise_positivity_times <- function(dt, ct_threshold = 37) {
  
  dt_proc <- data.table::copy(dt)
  
  # making sure we return the first time Ct values are closest to 37
  # (i.e. on the way up, not on the way down)
  # timing at which inidividuals begin to be detectable by PCR test
  pcr_pos_times_me_pre <- dt_proc[t < 6, .SD[which.min(abs(median - ct_threshold))],
                                  by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "median")]
  pcr_pos_times_lo_pre <- dt_proc[t < 6, .SD[which.min(abs(lo95 - ct_threshold))],
                                  by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "lo95")]
  pcr_pos_times_hi_pre <- dt_proc[t < 6, .SD[which.min(abs(hi95 - ct_threshold))],
                                  by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "hi95")]
  
  # merging into single data.table
  pcr_pos_times_pre <- merge(
    merge(pcr_pos_times_me_pre, pcr_pos_times_lo_pre,
          by = c("predictor", "regressor_category")),
    pcr_pos_times_hi_pre,
    by = c("predictor", "regressor_category"))
  
  # changing column names
  setnames(pcr_pos_times_pre, c("t", "t.x", "t.y"), c("t_hi95_pre", "t_me_pre", "t_lo95_pre"))
  
  # returning neat data.table
  pcr_pos_times_pre <- pcr_pos_times_pre[, c("predictor", "regressor_category", "t_me_pre", "t_lo95_pre", "t_hi95_pre")]
  
  # now making sure we return the second time Ct values are closest to 37
  pcr_pos_times_me_post <- dt_proc[t > 6, .SD[which.min(abs(median - ct_threshold))],
                                   by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "median")]
  pcr_pos_times_lo_post <- dt_proc[t > 6, .SD[which.min(abs(lo95 - ct_threshold))],
                                   by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "lo95")]
  pcr_pos_times_hi_post <- dt_proc[t > 6, .SD[which.min(abs(hi95 - ct_threshold))],
                                   by = c("predictor", "regressor_category")][, c("t", "predictor", "regressor_category", "hi95")]
  
  # merging
  pcr_pos_times_post <- merge(
    merge(pcr_pos_times_me_post, pcr_pos_times_lo_post,
          by = c("predictor", "regressor_category")),
    pcr_pos_times_hi_post,
    by = c("predictor", "regressor_category"))
  
  # changing column names
  setnames(pcr_pos_times_post, c("t", "t.x", "t.y"), c("t_lo95_post", "t_me_post", "t_hi95_post"))
  
  # returning neat data.table
  pcr_pos_times_post <- pcr_pos_times_post[, c("predictor", "regressor_category",
                                               "t_me_post", "t_lo95_post", "t_hi95_post")]
  
  out <- merge(pcr_pos_times_pre, pcr_pos_times_post, by = c("predictor", "regressor_category"))
  
  out <- out[!is.na(predictor)][order(regressor_category, predictor)]
  
  return(out)
  
}
