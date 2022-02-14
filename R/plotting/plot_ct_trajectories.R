# function for plotting the Ct trajectories. Needs the fitted Ct draws
# from the inference (ct_dt_draws) and the data.table of the data the 
# model was fit to (ct_dt)

plot_ct_trajectories <- function(ct_dt_draws, ct_dt) {
  
  p_out <- ct_dt_draws %>% 
    ggplot(aes(x = time)) +
    geom_line(aes(y = me)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_point(data = ct_dt, aes(x = t, y = ct_value, colour = factor(pcr_res))) + 
    facet_wrap(vars(factor(id))) +
    custom_plot_theme()
  
  return(p_out)
}
