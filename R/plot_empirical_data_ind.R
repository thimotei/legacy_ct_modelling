# quick plot of Ct data stratified by variant

plot_empiricial_data_ind <- function(dt_delta, dt_omicron) {
  
  p1 <- dt_delta %>%  
    ggplot(aes(x = t_first_test, y = ct_value)) +
    geom_point(aes(colour = factor(result))) +
    geom_line() +
    facet_wrap(~id) + 
    custom_plot_theme()
  
  p2 <- dt_omicron %>%  
    ggplot(aes(x = t_first_test, y = ct_value)) +
    geom_point(aes(colour = factor(result))) +
    geom_line() +
    facet_wrap(~id) + 
    custom_plot_theme()
  
  p_both <- p1 + p2 + plot_layout(guides = "collect")
  
  return(p_both)
}
