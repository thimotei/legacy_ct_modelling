extract_ct_trajectories <- function(fit, variable = "ct", inf_time = TRUE) {
  dt_draws <- fit$draws(variables = variable)
  dt_draws <- data.table::as.data.table(dt_draws)

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
