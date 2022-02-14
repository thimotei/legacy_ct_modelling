extract_ct_fits <- function(ct_dt_draws) {
  
  ct_dt_out <- ct_dt_draws[, c("id", "time") := tstrsplit(variable, ",")] %>% 
  .[, id := str_remove(id, "ct\\[")] %>% 
  .[, time := str_remove(time, "]")] %>% 
  .[, time := as.numeric(time)] %>% 
  .[, id := factor(id)] %>% 
  .[, c("time", "iteration", "chain", "id", "value")] %>% 
  .[order(id, time)] %>% 
  .[,  value_unscaled := (mx - mn) * value + mn]
  
  return(ct_dt_out)
}
