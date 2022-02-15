ct_trajectory_summarise <- function(ct_dt_draws, 
                                    by = c("id", "time")) {
  
  ct_dt_summary_out <- ct_dt_draws[, .(me = quantile(value_unscaled, c(0.5)),
                                       lo = quantile(value_unscaled, c(0.025)),
                                       hi = quantile(value_unscaled, c(0.975))), 
                                   by = by]
  return(ct_dt_summary_out)
}
