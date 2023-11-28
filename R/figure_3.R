# plotting functions specifically for figure 3

plot_ct_trajectory_panel <- function(pop_ct_draws_cens,
                                     pop_ct_draws_summ,
                                     regressor_category_arg,
                                     factor_arg,
                                     title_arg,
                                     no_draws = 1000,
                                     c_lod = 40,
                                     positivity_threshold = 37,
                                     lft_threshold = 30,
                                     alpha_arg = 0.2) {
  
  pop_ct_draws_cens <- pop_ct_draws_cens[
    regressor_category == regressor_category_arg][,
    predictor := forcats::fct_relevel(predictor, factor_arg)
    ]
  
  pop_ct_draws_summ <- pop_ct_draws_summ[
    regressor_category == regressor_category_arg][,
    predictor := forcats::fct_relevel(predictor, factor_arg)
    ]
  
  plot_out <- ggplot(data = pop_ct_draws_summ) + 
    geom_line(data = pop_ct_draws_cens[.draw <= no_draws],
              aes(x = t,
                  y = ct_value,
                  colour = predictor,
                  group = id),
              alpha = 0.002) +
    geom_line(aes(x = t,
                  y = median,
                  colour = predictor),
              linewidth = 0.3) +
    geom_ribbon(aes(x = t,
                    ymin = lo90,
                    ymax = hi90,
                    fill = predictor),
                alpha = alpha_arg) +
    geom_hline(aes(yintercept = c_lod), linetype = "dashed") +
    geom_text(aes(0,
                  c_lod,
                  label = "Limit of\ndetection", 
                  vjust = -0.2),
              colour = "black",
              size = 2.5) + 
    geom_hline(aes(yintercept = positivity_threshold),
               linetype = "dashed") +
    geom_text(aes(0,
                  positivity_threshold,
                  label = "PCR +ve\nthreshold", 
                  vjust = -0.2),
              size = 2.5) + 
    # very approximate lft threshold - check this with collaborators
    geom_hline(aes(yintercept = lft_threshold),
               linetype = "dashed") +
    geom_text(aes(0,
                  lft_threshold,
                  label = "LFT +ve\nthreshold",
                  vjust = -0.2),
              size = 2.5) +
    scale_y_reverse() + 
    coord_cartesian(clip = "off", ylim = c(40, NA)) + 
    labs(x = "Time since exposure (days)",
         y = "Ct value",
         title = title_arg) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    guides(fill = "none")
  
  return(plot_out)
  
}

plot_effect_panel <- function(param_pp,
                              regressor_category_arg,
                              baseline_arg,
                              factor_arg) {
  
  # creating some fake data to fix the limits on each facet, but changing
  # the limits between the different variables - bit hacky, but best way
  # I could find given the structure of the plotting code
  # blank_data <- data.table(
  #   variable = factor(rep(c("Ct value at peak", "Ct value at peak",
  #                           "Time of peak", "Time of peak",
  #                           "Time until LOD", "Time until LOD"),
  #                         13)),
  #   predictor = factor(c(
  #     rep("Delta", 6),
  #     rep("Omicron (BA.1)", 6),
  #     rep("Omicron (BA.2)", 6),
  #     rep("Symptomatic", 6),
  #     rep("Asymptomatic", 6),
  #     rep("3 exposures", 6),
  #     rep("4 exposures", 6),
  #     rep("5 exposures", 6),
  #     rep("6 exposures", 6),
  #     rep("Age: 20-34", 6),
  #     rep("Age: 35-49", 6),
  #     rep("Age: 50-64", 6),
  #     rep("Age: 65+", 6))),
  #   regressor_category = factor(c(
  #     rep("VOC", 18),
  #     rep("Symptom status", 12),
  #     rep("Number of exposures", 24),
  #     rep("Age", 24))),
  #   y = rep(c(14, 21, 4, 10, 15, 40), 13)
  # )
  
  # filtering by regressor category, to stratify within each panel
  param_pp_plot <- param_pp[
    regressor_category == regressor_category_arg][,
    predictor := forcats::fct_relevel(predictor, factor_arg)
    ]
  
  # blank_data_plot <- blank_data[regressor_category == regressor_category_arg]
  
  # the main plotting code
  plot_out <- ggplot() + 
    geom_pointrange(data = param_pp_plot,
                    aes(x = predictor,
                        y = median,
                        ymin = lo90,
                        ymax = hi90,
                        colour = predictor)) +
    geom_hline(data = param_pp_plot[predictor == baseline_arg],
               aes(yintercept = median), linetype = 2) +
    # geom_blank(data = blank_data_plot,
    #            aes(x = predictor, y = y)) +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    facet_rep_wrap(~variable,
                   scales = "free",
                   nrow = 3,
                   strip.position = "bottom") + 
    # scale_color_aaas() +
    custom_plot_theme() +
    labs(y = "Density") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  return(plot_out)
  
}

figure_3_subpanel <- function(p1, p2, rel_widths_arg) {
  
  plot_out <- plot_grid(p1, p2, rel_widths = rel_widths_arg)
  
  return(plot_out)
}

figure_3_t_lod <- function(dt_in, ct_threshold = 37) {
  
  dt_proc <- copy(dt_in)
  
  dt_out <- dt_proc[t > t_p, .SD[which.min(abs(ct_value - ct_threshold))],
                    by = c(".draw",
                           "predictor",
                           "regressor_category")]
  
  return(dt_out)
}

figure_3_t_lod_sum <- function(dt_in, ct_threshold = 37) {
  
  dt_proc <- copy(dt_in)
  
  dt_proc_cal <- dt_proc[t > t_p, .SD[which.min(abs(ct_value - ct_threshold))],
                         by = c(".draw",
                                "predictor",
                                "regressor_category")]
  
  dt_out <- dt_proc_cal[, .(variable = "Time until PCR-",
                            lo30 = quantile(t, 0.35), 
                            lo60 = quantile(t, 0.20),
                            lo90 = quantile(t, 0.05),
                            lo95 = quantile(t, 0.025),
                            median = quantile(t, 0.5),
                            hi30 = quantile(t, 0.65), 
                            hi60 = quantile(t, 0.8),
                            hi90 = quantile(t, 0.95),
                            hi95 = quantile(t, 0.975)),
                        by = c("predictor", "regressor_category")]
  
  return(dt_out)
}

add_baseline_effect_sizes <- function(dt_in) {
  # manually producing a single row with the same structure as the standard
  # effect size panel, to add the baseline draw
  dt_proc <- dt_in[, .(variable = "Time until PCR-",
                       lo30 = quantile(t, 0.35), 
                       lo60 = quantile(t, 0.20),
                       lo90 = quantile(t, 0.05),
                       lo95 = quantile(t, 0.025),
                       median = quantile(t, 0.5),
                       hi30 = quantile(t, 0.65), 
                       hi60 = quantile(t, 0.8),
                       hi90 = quantile(t, 0.95),
                       hi95 = quantile(t, 0.975)),
                   by = c("predictor", "regressor_category")]
  
  # same baseline description
  baseline_description <- "Baseline:
                         Omicron,
                         4 exposures,
                         symptomatic,
                         Age:35—49"
  
  # isolating the baseline category 
  dt_out <- dt_proc[
    predictor == "Omicron (BA.1)"
  ][, regressor_category := baseline_description]
  
  return(dt_out) 
}
