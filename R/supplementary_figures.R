time_between_last_neg_first_pos <- function(dt) {
  
  obs <- copy(dt)
  
  obs_diff <- rbind(
    obs[result == "Negative"][, .SD[swab_date == max(swab_date)], by = id], 
    obs[result == "Positive"][, .SD[swab_date == min(swab_date)], by = id])[
      order(id, result)
    ]
  
  obs_diff[, t_between_neg_pos := as.numeric(abs(swab_date - shift(swab_date)),
                                             units = "days"),
           by = id]
  
  p_out <- obs_diff[!is.na(t_between_neg_pos)] %>% 
    ggplot() + 
    geom_density(aes(x = t_between_neg_pos,
                     colour = t_between_neg_pos,
                     fill = t_between_neg_pos),
                 alpha = 0.3, fill = "#377EB8") + 
    custom_plot_theme() + 
    theme(legend.position = "none",
          panel.background = element_rect(colour = "black", fill = NA, size = 1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()) + 
    labs(x = "Time between last negative and first positive", y = "Density")
  
  return(p_out)
}

time_since_last_dose_plot <- function(dt) {
  
  obs <- data.table::copy(dt)
  
  obs[, iid := .GRP, c("id", "infection_id")][,
                                              no_vaccines := forcats::fct_relevel(no_vaccines, "2")]
  
  fig_2_spline_full <- obs[ct_value < 40] %>%
    ggplot(aes(x = time_since_last_dose*365, y = ct_value,
               colour = symptoms, fill = symptoms)) +
    geom_point(alpha = 0.5) +
    geom_smooth(alpha = 0.2) +
    # facet_rep_wrap(~) +
    scale_y_reverse() +
    theme(axis.line = element_line(),
          strip.background = element_rect(fill = "white", colour = "black", size = 1),
          strip.text = element_text(colour = "black"),
          legend.position = "none") +
    scale_colour_brewer(palette = "Set1",
                        labels = c("Symotomatic", "Asymptomatic")) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Symotomatic", "Asymptomatic")) +
    labs(x = "Time since last vaccine dose", y = "Ct value")
  
}

cumulative_shedding <- function(shedding_plot_dt,
                                ip_plot_dt,
                                regressor_category_arg,
                                factor_arg,
                                title_arg) {
  
  shedding_plot_dt <- shedding_plot_dt[regressor_category == regressor_category_arg]
  ip_plot_dt <- ip_plot_dt[regressor_category == regressor_category_arg]
  
  shedding_plot_dt <- shedding_plot_dt[,
                                       predictor := fct_relevel(predictor, factor_arg)
  ]
  
  ip_plot_dt <- ip_plot_dt[,
                           predictor := fct_relevel(predictor, factor_arg)
  ]
  
  p_out <- ggplot(data = shedding_plot_dt) +
    geom_line(aes(x = t, y = me, colour = predictor)) +
    geom_vline(aes(xintercept = ip_median, colour = predictor),
               linetype = "dashed") +
    geom_ribbon(aes(x = t, ymin = lo, ymax = hi, fill = predictor),
                alpha = 0.2) +
    geom_density(data = ip_plot_dt,
                 aes(x = inc_mean_nat, fill = predictor),
                 alpha = 0.1, colour = NA) +
    labs(x = "Time since exposure (days)",
         y = "Cumulative shedding (%)",
         title = title_arg) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    scale_y_continuous(labels = scales::percent) +
    lims(x = c(NA, 10))
  
  return(p_out)
  
}

shedding_effect_size_panel <- function(shedding_effects,
                                       regressor_category_arg,
                                       factor_arg,
                                       baseline_arg) {

  shedding_effects_plot <- shedding_effects[
    regressor_category == regressor_category_arg
  ]

  shedding_effects_plot <- shedding_effects_plot[,
  predictor := fct_relevel(predictor, factor_arg)
  ]

  shedding_effects_plot_baseline <- shedding_effects_plot[
    predictor == baseline_arg]

  p_out <- ggplot() +
    geom_pointrange(data = shedding_effects_plot,
                    aes(x = predictor,
                        y = me,
                        ymin = lo,
                        ymax = hi,
                        colour = predictor)) +
    geom_hline(data = shedding_effects_plot_baseline,
               aes(yintercept = me), linetype = 2) +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          legend.position = "none") +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0.025, 0.3)) +
    labs(y = "Cumulative shedding (%)")

  return(p_out)

}

cumulative_shedding_subpanel <- function(p1, p2, rel_widths_arg) {

  plot_both <- plot_grid(p1, p2, rel_widths = rel_widths_arg)

  plot_out <- plot_both
  
  return(plot_out)
}
