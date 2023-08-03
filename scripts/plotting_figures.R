library(cowplot)
library(ggridges)

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

pop_draws <- extract_ct_params(adj_draws, by = "predictor", mean = FALSE) 

pop_ct_draws <- pop_draws[.draw <= 1000] %>%
  transform_to_model() %>%
  simulate_cts(time_range = seq(0, 35, 0.25), obs_noise = FALSE)

pop_ct_draws_cens <- pop_ct_draws[ct_value <= 45]

sum_cols <- c("value", "t", by = "predictor")

pop_ct_sum <- summarise_draws(pop_ct_draws_cens[, value := ct_value][,
  ..sum_cols],
  by = setdiff(sum_cols, "value"))

p1 <- ggplot() + 
  geom_line(data = pop_ct_draws_cens, aes(x = t, y = ct_value, colour = predictor, group = id),
            alpha = 0.005) +
  geom_ribbon(data = pop_ct_sum, aes(x = t, ymin = lo90, ymax = hi90,
                  fill = predictor), alpha = 0.3) +
  geom_hline(data = pop_ct_sum, aes(yintercept = 40), linetype = "dashed") +
  facet_rep_wrap(vars(predictor), ncol = 1) +
  custom_plot_theme() +
  scale_y_reverse(limits = c(45, NA)) + 
  labs(x = "Time since inferred exposure (days)", y = "Ct value") + 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

param_pp <- pop_draws %>%
  transform_to_natural() %>% 
  melt_draws(ids = c(".chain", ".iteration", ".draw", by = "predictor")) %>% 
  .[!variable %in% c("t_s", "c_s")]

set_draw_limits <- function(draws, variable_arg, min, max) {
  draws[variable == variable_arg & 
    value > min & value < max]
}

c_p_draws <- set_draw_limits(param_pp, "c_p", min = 12 , max = 22)
t_p_draws <- set_draw_limits(param_pp, "t_p", min = 1, max = 6.5)
t_lod_draws <- set_draw_limits(param_pp, "t_lod", min = 18, max = 50)

param_pp <- rbind(c_p_draws, t_p_draws, t_lod_draws)

p2 <- param_pp[predictor == "Baseline (Omicron & symptomatic & 3 vaccines)",
               predictor := "Omicron"] %>% 
  update_variable_labels() %>%
  ggplot(aes(x = value, y = forcats::fct_rev(predictor),
             fill = predictor)) + 
  geom_density_ridges(scale = 4, alpha = 0.4) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(vars(variable), scales = "free", nrow = 3) +
  custom_plot_theme() + 
  labs(x = "Value", y = "")

cowplot::plot_grid(p1, p2)

