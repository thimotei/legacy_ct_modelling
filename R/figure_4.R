figure_4_data <- function(dt_in, ct_threshold) {
  
  dt_proc <- copy(dt_in)
  
  dt_proc <- dt_proc[t < t_p, .SD[which.min(abs(ct_value - ct_threshold))],
                    by = c(".draw",
                           "inc_mean_nat",
                           "predictor",
                           "regressor_category")]
  
  dt_out <- dt_proc[!is.na(predictor) == TRUE][predictor != "Asymptomatic"]
  
  return(dt_out)
}

figure_4_panel <- function(dt, regressor_category_arg,
                           title_arg,
                           factor_arg) {
  
  dt_proc <- data.table::copy(dt)
  
  dt_plot <- dt_proc[
    regressor_category == regressor_category_arg][, 
    predictor := forcats::fct_relevel(predictor, factor_arg)
  ]
  
  p_out <- dt_plot %>% 
    ggplot(aes(x = inc_mean_nat,
               y = t,
               colour = predictor)) + 
    geom_point(alpha = 0) + 
    geom_density_2d() + 
    theme_minimal() +
    labs(x = "Incubation period mean posterior draw",
         y = "Time threshold reached (Ct = 20)",
         title = title_arg) + 
    theme(legend.position = "bottom",
          legend.title = element_blank()) + 
    lims(x = c(3, 8),
         y = c(3, 8)) 
  
  return(p_out)
}

add_marginals <- function(p_plot) {
  
  p_out <- ggMarginal(p_plot,
                      groupColour = TRUE,
                      groupFill = TRUE,
                      alpha = 0.1)
  
  return(p_out)
}
