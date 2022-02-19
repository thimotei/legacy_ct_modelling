extract_posterior_predictions <- function(fit, obs) {
  dt_draws <- fit$draws(variables = "sim_ct")
  dt_draws <- data.table::as.data.table(dt_draws)

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
