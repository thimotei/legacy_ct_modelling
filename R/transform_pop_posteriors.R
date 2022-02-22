transform_pop_posteriors <- function(pop_params, draws) {
  
  draws_long <- draws[variable %in% pop_params] %>%
    .[, chain := NULL] %>% 
    .[, iteration := as.numeric(iteration)] %>% 
    dcast(., iteration + voc ~ variable, 
          value.var = "value", 
          fun.aggregate = function(x) x[1]) %>% 
    .[order(voc, iteration)]
  
  draws_wide_abs <- draws_long[, c_0_abs := c_0,
                               by = c("iteration", "voc")] %>%
    .[, c_s_mean_abs := c_0*plogis(c_s_mean), 
      by = c("iteration", "voc")] %>% 
    .[, c_p_mean_abs := plogis(c_p_mean)*c_s_mean_abs, 
      by = c("iteration", "voc")] %>% 
    .[, t_p_mean := exp(t_p_mean), 
      by = c("iteration", "voc")] %>% 
    .[, t_s_mean:= exp(t_s_mean), 
      by = c("iteration", "voc")] %>% 
    .[, t_lod_mean := exp(t_lod_mean), 
      by = c("iteration", "voc")] %>% 
    .[, t_p_mean_abs := t_p_mean, 
      by = c("iteration", "voc")] %>% 
    .[, t_s_mean_abs := t_p_mean_abs + t_s_mean, 
      by = c("iteration", "voc")] %>% 
    .[, t_lod_mean_abs :=  t_p_mean + t_s_mean + t_lod_mean, 
      by = c("iteration", "voc")] %>% 
    .[order(voc, iteration)]
  
  return(draws_wide_abs)
}
