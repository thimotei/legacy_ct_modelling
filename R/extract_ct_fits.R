extract_ct_fits <- function(dt_draws) {

  ct_dt_out <- dt_draws[variable %like% "ct"][,
   c("id", "time") := tstrsplit(variable, ",")
  ][,
   id := str_remove(id, "ct\\[")][,
   time := str_remove(time, "\\]")][,
   time := as.numeric(time)][,
   id := factor(id)][,
   c("time", "iteration", "chain", "id", "value")][
   order(id, time)]

  inf_time_draws <- dt_draws[variable  %like% "T_e"][,
    id := str_remove(variable, "T_e\\[")][,
    id := str_remove(id, "\\]")][,
    .(id, inf_time = value, iteration, chain)]


  ct_dt_out <- ct_dt_out[inf_time_draws, on = c("id", "iteration", "chain")]
  ct_dt_out[, time_since_first_pos := time - inf_time]
  return(ct_dt_out[])
}
