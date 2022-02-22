plot_pop_posteriors <- function(input_data) {
  
  p_out <- input_data %>% 
    ggplot() + 
    geom_density(aes(x = value, fill = VOC), alpha = 0.2) +
    facet_wrap(~variable, nrow = 2, scales = "free") + 
    labs(x = "Value", y = "Probability density")
  
  return(p_out)
}
