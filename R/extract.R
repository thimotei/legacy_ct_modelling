extract_draws <- function(fit, params = c("c_0", "c_p_mean",
                                          "c_s_mean", "t_p_mean",
                                          "t_s_mean", "t_lod_mean")) {
  draws <- fit$draws(format = "df", variables = params)
  draws <- data.table::as.data.table(draws)
  return(draws[])
}

melt_draws <- function(draws, ids = c(".chain", ".iteration", ".draw")) {
  data.table::melt(draws, id.vars = ids)
}

extract_ct_trajectories <- function(fit, variable = "ct", inf_time = TRUE) {
  dt_draws <- extract_draws(fit, params = variable)

  ct_dt_out <- dt_draws[,
   c("id", "time") := tstrsplit(variable, ",")
  ][,
   id := str_remove(id, paste0("ct", "\\["))][,
   time := str_remove(time, "\\]")][,
   time := as.numeric(time)][,
   id := factor(id)][,
   c("time", "iteration", "chain", "id", "value")][
   order(id, time)]

  if (inf_time) {
    inf_time_draws <- data.table::as.data.table(
      fit$draws(variables = "T_e")
    )[,
      id := str_remove(variable, "T_e\\[")][,
      id := str_remove(id, "\\]")][,
      .(id, inf_time = value, iteration, chain)]


    ct_dt_out <- ct_dt_out[inf_time_draws, on = c("id", "iteration", "chain")]
  }
  ct_dt_out[, time_since_first_pos := time - inf_time]
  return(ct_dt_out[])
}

extract_posterior_predictions <- function(fit, obs) {
  dt_draws <- extract_draws(fit, "sim_ct")

  simulated_cts <- dt_draws[,
    obs := str_remove(variable, "sim_ct\\[")][,
    obs := str_remove(obs, "\\]")][,
    .(obs = as.numeric(obs), sim_ct = value, iteration = as.numeric(iteration),
      chain = as.numeric(chain))
  ][order(obs)]

  if (!missing(obs)) {
    simulated_cts <- merge(
      obs[, obs := 1:.N], simulated_cts, by = "obs"
    )
  }
  return(simulated_cts[])
}
