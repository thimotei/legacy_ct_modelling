extract_subjects <- function(dt) {
  subjects <- data.table::copy(dt)
  subjects <- subjects[, .SD[1,], by = "id"]
  
  return(subjects[])
}

extract_draws <- function(fit, params = NULL, format = "df") {
  draws <- fit$draws(format = format, variables = params)
  draws <- data.table::as.data.table(draws)
  
  return(draws[])
}

extract_params <- function(draws, params, by) {
  if (!missing(by)) {
    params <- c(params, by)
  }
  cols <- intersect(colnames(draws), params)
  cols <- c(".iteration", ".draw", ".chain", cols)
  draws <- draws[, ..cols]
  
  return(draws[])
}

extract_ct_params <- function(draws, params = c("c_0", "c_p_mean",
                                                "c_s_mean[1]", "t_p_mean",
                                                "t_s_mean[1]", "t_lod_mean"),
                              mean = TRUE, by) {
  if (!mean) {
    params <- stringr::str_remove(params, "_mean")
    params <- stringr::str_remove(params, "\\[1\\]")
  }
  draws <- extract_params(draws, params = params, by)
  colnames(draws) <- stringr::str_remove(colnames(draws), "_mean")
  colnames(draws) <- stringr::str_remove(colnames(draws), "\\[1\\]")
  
  return(draws[])
}

extract_ip_params <- function(draws, params = c("inc_mean[1]", "inc_mean",
                                                "inc_sd[1]", "inc_sd"),
                              by) {
  draws <- extract_params(draws, params = params, by)
  colnames(draws) <- purrr::map_chr(
    colnames(draws), ~ stringr::str_split(., "\\[[0-9]\\]")[[1]][1]
  )
  
  return(draws)
}

extract_param_draws <- function(draws, onsets = TRUE) {
  
  if(onsets == TRUE) {
    draws <- cbind(
      extract_ct_params(draws),
      extract_ip_params(draws)[, .(inc_mean, inc_sd)]
    )
  }
  else {
    draws <- extract_ct_params(draws)
  }
  
  return(draws)
}

extract_coeffs <- function(draws, exponentiate = FALSE, design, variables) {
  beta_cols <- grep("beta_", colnames(draws), value = TRUE)
  cols <- c(".iteration", ".draw", ".chain", beta_cols)
  draws <- draws[, ..cols]
  
  if (exponentiate) {
    draws[, (beta_cols) := lapply(.SD, exp), .SDcols = beta_cols]
  }
  
  setnames(draws, beta_cols, gsub("beta_", "", beta_cols))
  
  draws <- melt_draws(draws)
  
  draws <- draws[,
                 coeff := as.numeric(stringr::str_extract(variable, "[0-9]+"))
  ][,
    variable := purrr::map_chr(
      variable, ~ stringr::str_split(., "\\[[0-9]+\\]")[[1]][1]
    )
  ]
  
  if (!missing(variables)) {
    draws <- draws[variable %in% variables]
  }
  
  if (!missing(design)) {
    design <- data.table::data.table(predictor = colnames(design))
    design <- design[!predictor %in% "(Intercept)"]
    design[, coeff := 1:.N]
    draws <- draws[design, on = "coeff"]
  }
  return(draws[])
}

melt_draws <- function(draws, ids = c(".chain", ".iteration", ".draw")) {
  data.table::melt(draws, id.vars = ids)
}

extract_ct_trajectories <- function(fit, variable = "ct", inf_time = TRUE) {
  dt_draws <- extract_draws(fit, params = variable, format = "array")
  
  obs_out <- dt_draws[, c("id", "time") := tstrsplit(variable, ",")
  ][, id := str_remove(id, paste0("ct", "\\["))
  ][, time := str_remove(time, "\\]")
  ][, time := as.numeric(time)
  ][, id := factor(id)
  ][, c("time", "iteration", "chain", "id", "value")
  ][order(id, time)]
  
  if (inf_time) {
    inf_time_draws <- extract_draws(
      fit, params = "t_inf", format = "array"
      )[, id := str_remove(variable, "t_inf\\[")
      ][, id := str_remove(id, "\\]")
      ][,.(id, inf_time = value, iteration, chain)]
    
    obs_out <- obs_out[inf_time_draws, on = c("id", "iteration", "chain")]
  }
  obs_out[, time_since_first_pos := time - inf_time]
  cols <- c("chain", "iteration")
  obs_out[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  return(obs_out[])
}

extract_pop_ct_trajectories <- function(fit,
                                        no_draws = 1000,
                                        tmin = 0,
                                        tmax = 30,
                                        tstep = 0.1,
                                        lower_ct_limit = 40,
                                        separate_baseline_covariates = TRUE,
                                        baseline_tag = "Baseline",
                                        other_covariates = TRUE,
                                        censor_output = TRUE,
                                        onsets_flag = onsets_flag) {
  # Extract posterior predictions
  draws <- extract_draws(fit)
  
  adj_draws <- adjust_params(draws, 
                             design = ct_model$design, 
                             onsets_flag = onsets_flag) 
  
  adj_draws <- adj_draws %>% update_predictor_labels()
  
  # adding baseline draws
  if(separate_baseline_covariates == TRUE & other_covariates == TRUE) {
    adj_draws <- add_baseline_to_draws(
      adj_draws, "Omicron (BA.1)", onsets_flag = onsets_flag)
    
    adj_draws <- add_baseline_to_draws(
      adj_draws, "4 exposures", onsets_flag = onsets_flag)
    
    adj_draws <- add_baseline_to_draws(
      adj_draws, "Symptomatic", onsets_flag = onsets_flag)
    
    adj_draws <- add_baseline_to_draws(
      adj_draws, "Age: 35-49", onsets_flag = onsets_flag)
    
  } else if(separate_baseline_covariates == FALSE & other_covariates == FALSE) {
    adj_draws[is.na(predictor), predictor := "Omicron (BA.1)"]
  } 
  
  # simulating Ct trajectories
  pop_ct_draws <- adj_draws %>%
    transform_to_model(., onsets_flag = onsets_flag) %>%
    simulate_cts(time_range = seq(tmin, tmax, tstep), 
                 obs_noise = FALSE)
  
  # returning the number of draws set in function call
  pop_ct_draws <- pop_ct_draws[, .SD[.draw %in% 1:no_draws],
                               by = "predictor"]
  
  # adding regressor categories
  pop_ct_draws <- add_regressor_categories(pop_ct_draws)
  
  # censoring outputs to the lower limit of detection (Ct = 40) for
  # simpler interpretability
  out <- pop_ct_draws[ct_value > lower_ct_limit, ct_value := 40]
  
  return(out)
  
}

extract_posterior_predictions <- function(fit, obs) {
  dt_draws <- extract_draws(fit, "sim_ct", format = "array")
  
  simulated_cts <- dt_draws[, obs := str_remove(variable, "sim_ct\\[")
                            ][, obs := str_remove(obs, "\\]")
                            ][, .(obs = as.numeric(obs), 
                                  sim_ct = value, 
                                  iteration = as.numeric(iteration),
                                  chain = as.numeric(chain))
                            ][order(obs)]
  
  if (!missing(obs)) {
    simulated_cts <- merge(
      obs[order(id), obs := 1:.N], simulated_cts, by = "obs"
    )
  }
  return(simulated_cts[])
}

extract_shedding <- function(dt,
                             vl_flag = FALSE,
                             trim_flag = TRUE,
                             pcr_pos_threshold = 37) {
  
  out <- data.table::copy(dt)
  
  out <- out[ct_value < pcr_pos_threshold]
  
  # calculate the cumulative area under Ct curves, as proportion of total
  # area, by ID (keeping predictor labels too)
  if(vl_flag == TRUE) {
    
    # converting to viral load from Ct value 
    out[, vl := 2^(pcr_pos_threshold - ct_value),
        by = c(".draw", "predictor", "id")]
    
    # calculating total area under curve for each trajectory
    total_auc_dt <- out[ct_value <= pcr_pos_threshold,
                        .(total_auc = pracma::trapz(t, vl)),
                        c(".draw", "predictor", "id")]
    
    # merging with viral load trajectories
    out <- merge(out,
                 total_auc_dt,
                 by = c(".draw", "predictor", "id"))
    
    # calculating AUC under VL curve
    out[ct_value <= pcr_pos_threshold,
        auc_prop := pracma::cumtrapz(t, vl)/pracma::trapz(t, vl),
        by = c(".draw", "predictor", "id")]
  } else {
    
    # calculating total area under curve for each trajectory
    total_auc_dt <- out[, .(total_auc = pracma::trapz(t, ct_value)),
                        c(".draw", "predictor", "id")]
    
    # merging with viral load trajectories
    out <- merge(out,
                 total_auc_dt,
                 by = c(".draw", "predictor", "id"))
    
    # calculating AUC under Ct curve
    out[ct_value <= pcr_pos_threshold,
        auc_prop := pracma::cumtrapz(t, ct_value)/total_auc,
        by = c(".draw", "predictor", "id")]
  }
  
  # adding the categories for the many predictors
  add_regressor_categories(out)
  
  if(trim_flag == TRUE) {
    out <- out[, c("t", "ct_value", "vl", "auc_prop",
                   ".draw", "id", "predictor",
                   "regressor_category")]
  }
  
  # shift time back a day, to account for the onset times being shifted
  # back one day
  # out <- out[, t := t - 1]
  
  return(out)
}

extract_ip_draws <- function(draws,
                             factor_order,
                             keep_asymptomatic = FALSE,
                             trim_flag = TRUE) {
  
  # extracting the adjusted draws
  adj_draws <- adjust_params(draws,
                             design = ct_model$design,
                             onsets_flag = TRUE) %>%
    update_predictor_labels()
  
  adj_draws <- add_baseline_to_draws(adj_draws, "Omicron (BA.1)", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "4 exposures", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "Symptomatic", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "Age: 35-49", onsets_flag = TRUE)
  
  # add_baseline_to_draws(cumulative_shedding, "Omicron (BA.1)", onsets_flag = TRUE)
  
  # calculating draws from the incubation period using the adjusted draws
  out <- extract_ip_params(adj_draws, by = "predictor")[,
                                                        ip_draw := rlnorm(1, inc_mean, inc_sd), by = c(".draw", "predictor")
  ][, inc_mean_nat := exp(inc_mean)]
  
  # remove the few negative draws (where are these coming from?)
  out <- out[!is.na(ip_draw)]
  
  out[,  predictor := fct_relevel(predictor,
                                  factor_order)
  ]
  
  if(keep_asymptomatic == FALSE) {
    # removing artificial incubation period estimate for asymptomatic individuals
    out <- out[predictor != "Asymptomatic"]
  }
  
  # adding regressor categories
  add_regressor_categories(out)
  
  if(trim_flag == TRUE) {
    out <- out[, c(".draw", "inc_mean", "inc_mean_nat", "inc_sd", 
                   "ip_draw", "predictor", "regressor_category")]
  }
  
  return(out)
  
}