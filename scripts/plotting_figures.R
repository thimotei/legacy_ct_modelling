library(cowplot)
library(ggridges)
library(lemon)

# Extract posterior predictions
draws <- extract_draws(fit)

adj_draws <- adjust_params(draws, design = ct_model$design, inc_period = TRUE) %>%
  update_predictor_labels()

# Filter for just adjustments that summary shows appear to differ from base case
# adj_draws <- adj_draws[
#   predictor %in% c("no_vaccines2", "symptomsasymptomatic", 
#                    "VOCDelta", "VOCBA.2", "time_since_last_dose")
# ] 

adj_draws <- rbind(
  extract_param_draws(draws, inc_period = TRUE)[,
  predictor := "Omicron (baseline)"
  ],
  adj_draws
)[, predictor := factor(predictor)]

custom_plot_theme <- function() {
  
  custom_plot_theme = list(
    theme_cowplot(11),
    theme(axis.line = element_line()),
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.position = "none"))
  
  return(custom_plot_theme)
}

#--- figure 2 (maybe figure 2a)

make_figure_2 <- function(adj_draws) {
  
  pop_draws <- extract_ct_params(adj_draws, by = "predictor", mean = FALSE) 
  
  pop_ct_draws <- pop_draws %>%
    transform_to_model() %>%
    simulate_cts(time_range = seq(0, 15, 0.25), obs_noise = FALSE)
  
  pop_ct_draws_cens <- pop_ct_draws[ct_value <= 45]
  
  sum_cols <- c("value", "t", by = "predictor")
  
  pop_ct_sum <- summarise_draws(pop_ct_draws_cens[, value := ct_value][,
                                                                       ..sum_cols],
                                by = setdiff(sum_cols, "value"))
  
  pop_ct_draws_plot <- pop_ct_draws_cens[, .SD[.draw %in% 1:500],
                                         by = "predictor"]
  p21 <- ggplot() + 
    geom_line(data = pop_ct_draws_plot,
              aes(x = t, y = ct_value, colour = predictor, group = id),
              alpha = 0.002) +
    geom_ribbon(data = pop_ct_sum, aes(x = t, ymin = lo90, ymax = hi90,
                                       fill = predictor), alpha = 0.2) +
    geom_hline(data = pop_ct_sum, aes(yintercept = 40), linetype = "dashed") +
    facet_rep_wrap(vars(predictor), ncol = 1) +
    custom_plot_theme() +
    scale_y_reverse(limits = c(45, NA)) + 
    labs(x = "Time since exposure (days)", y = "Ct value")
  # scale_colour_brewer(palette = "Accent")
  
  param_pp <- pop_draws %>%
    transform_to_natural() %>% 
    melt_draws(ids = c(".chain", ".iteration", ".draw", by = "predictor")) %>% 
    .[!variable %in% c("t_s", "c_s")]
  
  set_draw_limits <- function(draws, variable_arg, min, max) {
    draws[variable == variable_arg & 
            value > min & value < max]
  }
  
  c_p_draws <- set_draw_limits(param_pp, "c_p", min = 13.5 , max = 20.5)
  t_p_draws <- set_draw_limits(param_pp, "t_p", min = 1, max = 6.5)
  t_lod_draws <- set_draw_limits(param_pp, "t_lod", min = 30, max = 55)
  
  param_pp <- rbind(c_p_draws, t_p_draws, t_lod_draws)
  
  p22 <- param_pp %>% 
    update_variable_labels() %>%
    ggplot(aes(x = value, y = forcats::fct_rev(predictor),
               fill = predictor)) + 
    geom_density_ridges(scale = 4, alpha = 0.4) +
    theme(legend.position = "none") +
    # scale_fill_brewer(palette = "Accent") +
    facet_wrap(vars(variable), scales = "free", nrow = 3) +
    custom_plot_theme() + 
    labs(x = "Value", y = "")
  
  p_2a <- cowplot::plot_grid(p21, p22)
  
  # ggsave(
  #   "outputs/figures/figure_2.png",
  #   p_all, width = 8, height = 12, bg = "white"
  # )
  
  #--- figure 3 (maybe figure 2b)
  # plotting posterior predictive
  
  ip_draws <- extract_ip_params(adj_draws, by = "predictor")
  
  pop_ip_draws <- ip_draws %>%
    simulate_ips(time_range = seq(0, 15, 0.05))
  
  sum_cols <- c("value", "t", by = "predictor")
  
  pop_ip_sum <- summarise_draws(
    pop_ip_draws[, ..sum_cols],
    by = setdiff(sum_cols, "value")
  )
  
  p31 <- ggplot() + 
    # geom_line(data = pop_ip_draws[, .SD[.draw %in% 1:500], by = "predictor"],
    #           aes(x = t, y = value, colour = predictor, group = id),
    #           alpha = 0.002) +
    geom_ribbon(data = pop_ip_sum, aes(x = t, ymin = lo90, ymax = hi90,
                                       fill = predictor), alpha = 0.5) +
    facet_rep_wrap(vars(predictor), ncol = 1) +
    custom_plot_theme() +
    lims(x = c(0, 9), y = c(0, 1)) + 
    labs(x = "Time since exposure (days)", y = "Probability of symptom onset")
  
  ip_draws_natural <- ip_draws %>% 
    transform_ip_to_natural() %>%
    melt(., id.vars = c(".iteration", ".draw", ".chain", "predictor"),
         measure.vars = c("nat_inc_mean", "nat_inc_sd"))
  
  # plotting parameter posteriors
  inc_mean_draws <- set_draw_limits(ip_draws_natural,
                                    "nat_inc_mean", min = 3 , max = 7)
  
  inc_sd_draws <- set_draw_limits(ip_draws_natural,
                                  "nat_inc_sd", min = 0, max = 4)
  
  ip_param_pp <- rbind(inc_mean_draws, inc_sd_draws)
  
  p32 <- ip_param_pp %>% 
    update_variable_labels() %>%
    ggplot(aes(x = value, y = forcats::fct_rev(predictor),
               fill = predictor)) + 
    geom_density_ridges(scale = 4, alpha = 0.4) +
    theme(legend.position = "none") +
    # scale_fill_brewer(palette = "Accent") +
    facet_wrap(vars(variable), scales = "free", nrow = 3) +
    custom_plot_theme() + 
    labs(x = "Value", y = "")
  
  p_2b <- cowplot::plot_grid(p31, p32)
  
  p_figure_2 <- cowplot::plot_grid(p_2a, p_2b)
  
  return(p_figure_2)
  
}

figure_2 <- make_figure_2(adj_draws)

ggsave("outputs/figures/figure_2.png",
       figure_2,
       bg = "white",
       width = 14,
       height = 11)

#--- figure 3 - difference between timing of peak and incubation period
# between at least delta and omicron
t_p_to_inc <- merge(param_pp[variable == "t_p"],
                    ip_draws,
                    by = c(".iteration", ".draw", ".chain", "predictor")) 

t_p_to_inc_sample <- t_p_to_inc[,
  inc_draw := rlnorm(1, inc_mean, inc_sd),
  by = c(".iteration", ".draw", "predictor")][,
  diff := inc_draw - value,
  by = c(".iteration", ".draw", "predictor")]

t_p_to_inc_sample %>% 
  ggplot() +
  geom_density(aes(x = diff, fill = predictor)) + 
  facet_wrap(vars(predictor)) +
  theme(legend.position = "none") +
  # scale_fill_brewer(palette = "Accent") +
  facet_rep_wrap(vars(predictor), nrow = 5) +
  custom_plot_theme() + 
  labs(x = "Value", y = "Difference between incubation period and peak timing") +
  lims(x = c(-1, 3))

t_p_to_inc_sample

#--- probability of detection results
pop_ct_draws <- pop_draws %>%
  transform_to_model() %>%
  simulate_cts(time_range = seq(0, 40, 0.25), obs_noise = FALSE)

pop_ct_draws[
  ct_value < 40, detect := 1][
  ct_value >= 40, detect := 0]

not_detected_sum <- pop_ct_draws[detect == 0, .(not_detected = .N), by = c("t", "predictor")]
detected_sum <- pop_ct_draws[detect == 1, .(detected = .N), by = c("t", "predictor")]
total <- pop_ct_draws[, .(total = .N), by = c("t", ".draw", "predictor")]

prob_detection_dt <- merge(detected_sum, not_detected_sum, all.x = TRUE)[,
  prob_detection := detected/4000]

prob_detecion_plot <- prob_detection_dt %>% 
  ggplot() + 
  geom_line(aes(x = t, y = prob_detection, colour = predictor)) +
  custom_plot_theme() + 
  labs(x = "Time (days) since exposure", y = "P(detection | infected)") 

virus_shed <- cbind(prob_detection_dt[, .(t), by = "predictor"][, .(t)], 
      prob_detection_dt[, .(shedding = 1 - cumsum(prob_detection)/sum(prob_detection)),
                        by = c("predictor")]) %>% 
  ggplot() + 
  geom_line(aes(x = t, y = shedding, colour = predictor)) +
  custom_plot_theme() + 
  labs(x = "Time (days) since exposure", y = "Viral shedding") +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

panel_r <- cowplot::plot_grid(prob_detecion_plot, virus_shed, nrow = 2)
panel_both <- cowplot::plot_grid(t_p_to_inc_plot, panel_r, 
                                 rel_widths = c(1, 1.25))

eff_plot <- plot_effect_summary(
  draws, ct_design = ct_model$design, variables = adj_params,
  adjustment_design = adjustment_model$design,
  variable_labels = update_variable_labels,
  col = predictor, position = position_dodge(width = 0.6)
) &
  scale_colour_brewer(palette = "Dark2") &
  labs(col = "Adjustment")

plot_effect_summary <- function(draws, ct_design, adjustment_design, variables,
                                variable_labels = function(dt, ...) {
                                  return(dt)
                                }, ...) {
  ct_plot <- draws %>%
    summarise_effects(design = ct_design, variables = variables) %>%
    update_predictor_labels() %>%
    variable_labels(reverse = TRUE) %>%
    plot_effects(...)
  
  adjustment_plot <- draws %>%
    summarise_adjustment(
      design = adjustment_design
    ) %>%
    update_predictor_labels() %>%
    variable_labels(reverse = TRUE) %>%
    plot_effects(
      trans = "identity", ...
    )
  
  summary <- ct_plot 
  
  return(summary)
}


figure_3 <- cowplot::plot_grid(eff_plot, panel_both)

ggsave("outputs/figures/figure_3.png", figure_3,
       bg = "white", width = 14, height = 8)

# t_p_to_inc_plot + prob_detecion_plot/virus_shed + plot_layout(guides = "collect")
