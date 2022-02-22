pop_pp_samples <- function(pop_posteriors_wide, t_max, t_step) {
  
  pop_pp_samples_out <- pop_posteriors_wide %>% 
  .[rep(.[, .I], (t_max)/t_step + 1)] %>%
  .[order(iteration, voc)] %>% 
  .[, time := seq(0, t_max, t_step), by = c("iteration", "voc")] %>% 
  .[, value := ct_hinge_long(time, c0 = c_0_abs, cp = c_p_mean_abs,
                             cs = c_s_mean_abs, clod = c_0_abs, 
                             te = 0, tp = t_p_mean, ts = t_s_mean, 
                             tlod = t_lod_mean),
    by = c("time", "iteration", "voc")]
  
  return(pop_pp_samples_out)
}
