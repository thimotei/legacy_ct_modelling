munge_figure_4_data <- function(draws,
                                dt_pop_ct_draws,
                                ct_threshold_arg) {
  
  adj_draws <- adjust_params(draws,
                             design = ct_model$design,
                             onsets_flag = TRUE) %>%
    update_predictor_labels()
  
  adj_draws <- add_baseline_to_draws(adj_draws, "Omicron (BA.1)", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "4 exposures", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "Symptomatic", onsets_flag = TRUE)
  adj_draws <- add_baseline_to_draws(adj_draws, "Age: 35-49", onsets_flag = TRUE)
  
  adj_draws_natural <- adj_draws %>%
    transform_to_model() %>%
    add_regressor_categories() %>% 
    .[!is.na(predictor) & predictor != "Asymptomatic"] %>% 
    .[, inc_mean_nat := exp(inc_mean), by = .draw]
  
  dt_ct_threshold <- calculate_ct_threshold(dt_pop_ct_draws,
                                            ct_threshold = ct_threshold_arg,
                                            trim_flag = TRUE)
  
  out <- merge(adj_draws_natural,
               dt_ct_threshold,
        by = c(".draw",
               "predictor",
               "regressor_category"))
  
  
  return(out)

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
    lims(x = c(3, 7),
         y = c(3, 7)) 
  
  return(p_out)
}

add_marginals <- function(p_plot) {
  
  p_out <- ggMarginal(p_plot,
                      groupColour = TRUE,
                      groupFill = TRUE,
                      alpha = 0.1)
  
  return(p_out)
}
