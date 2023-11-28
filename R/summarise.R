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
    ct_pp[, value := sim_ct], by = c("id", 
                                     "t_first_test_since_inf",
                                     "obs",
                                     "VOC",
                                     "ct_type")
  )
  return(ct_pp)
}

summarise_draws <- function(draws, by = c("id", "time")) {

  out <- draws[,
      .(lo30 = quantile(value, c(0.35), na.rm = TRUE),
        lo60 = quantile(value, c(0.20), na.rm = TRUE),
        lo90 = quantile(value, c(0.10), na.rm = TRUE),
        lo95 = quantile(value, c(0.05), na.rm = TRUE),
        median = quantile(value, c(0.5), na.rm = TRUE),
        hi30 = quantile(value, c(0.65), na.rm = TRUE),
        hi60 = quantile(value, c(0.80), na.rm = TRUE),
        hi90 = quantile(value, c(0.90), na.rm = TRUE),
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

summarise_effect_sizes_natural <- function(draws,
                                           ct_model = ct_model,
                                           add_baseline_flag = TRUE,
                                           onsets_flag = TRUE) {
  adj_params <- c("c_p", "t_p", "t_lod")
  
  if(add_baseline_flag == TRUE) {
    
    baseline_predictors <- c("Omicron",
                             "Symptomatic",
                             "4 exposures",
                             "Age: 35-49")
    
    adj_draws <- adjust_params(draws, 
                               design = ct_model$design,
                               onsets = onsets_flag) 
    
    adj_draws <- adj_draws %>% 
      update_predictor_labels()
    
    adj_draws <- add_baseline_to_draws(
      adj_draws, "Omicron (BA.1)", onsets_flag = onsets_flag)
  
    baseline_description <- "Baseline:
                             Omicron,  
                             4 exposures,
                             symptomatic,
                             Age:35—49"
    
    # adding baseline case with long description, for legend and vertical lines
    adj_draws <- add_baseline_to_draws(adj_draws,
                                       baseline_description,
                                       onsets_flag = onsets_flag)
    
    # adding the baseline results again for each subcategory, to plot
    # the results separately for each regressor category
    adj_draws <- add_baseline_to_draws(adj_draws, "4 exposures", onsets_flag = onsets_flag)
    adj_draws <- add_baseline_to_draws(adj_draws, "Symptomatic", onsets_flag = onsets_flag)
    adj_draws <- add_baseline_to_draws(adj_draws, "Age: 35-49", onsets_flag = onsets_flag)
  } else {
    
    baseline_predictors <- c("Omicron")
    
    adj_draws <- adjust_params(draws, 
                               design = ct_model$design, 
                               onsets = onsets_flag) 
    
    adj_draws <- adj_draws %>% 
      update_predictor_labels()
    
    adj_draws <- add_baseline_to_draws(
      adj_draws,
      "Omicron (BA.1)", 
      onsets_flag = onsets_flag)
  }
  
  pop_draws <- extract_ct_params(adj_draws, by = "predictor", mean = FALSE) %>% 
    transform_to_natural()
  
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

summarise_ct_threshold <- function(dt_in) {
  
  dt_out <- dt_in[, .(t_me = quantile(t, 0.5),
                      t_lo = quantile(t, 0.0025),
                      t_hi = quantile(t, 0.975)),
                  by = c("predictor",#
                         "regressor_category",
                         "direction")][
  order(direction, regressor_category, predictor)]
  
  return(dt_out)
}

summarise_posterior_differences <- function(adj_draws, param_arg) {
  
  by_vars = c(param_arg, ".draw", "predictor")
  
  dt_pop_long <- adj_draws[, ..by_vars]
  dt_pop_wide <- dcast(dt_pop_long,
                       .draw ~ predictor, 
                       value.var = param_arg)
  
  # VOC posterior differences
  dt_pop_diff_wide <- dt_pop_wide[,
  `:=` (diff_baseline_delta = baseline - Delta,
        diff_baseline_ba2 = baseline - `Omicron (BA.2)`,
        diff_baseline_asymp = baseline - Asymptomatic,
        diff_baseline_3_exps = baseline - `3 exposures`,
        diff_baseline_5_exps = baseline - `5+ exposures`,
        diff_baseline_20_34 = baseline - `Age: 20-34`,
        diff_baseline_50 = baseline - `Age: 50+`,
        diff_3_exps_5_exps = `3 exposures` - `5+ exposures`),
                                  by = ".draw"]
  
  dt_pop_diff_long <- melt(dt_pop_diff_wide[, c(".draw",
                                                "diff_baseline_delta",
                                                "diff_baseline_ba2",
                                                "diff_baseline_asymp",
                                                "diff_baseline_3_exps",
                                                "diff_baseline_5_exps",
                                                "diff_baseline_20_34",
                                                "diff_baseline_50",
                                                "diff_3_exps_5_exps")],
                           id.vars = ".draw")
  
  dt_out <- dt_pop_diff_long[, .(me = signif(quantile(value, 0.5), 2),
                                 lo = signif(quantile(value, 0.025), 2),
                                 hi = signif(quantile(value, 0.975), 2)),
                             by = "variable"]
  
  return(dt_out)
  
}
