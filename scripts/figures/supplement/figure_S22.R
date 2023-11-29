library(ggridges)
library(gridExtra)

# Sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# Either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
fit_main <- readRDS("outputs/fits/main.rds")
draws <- extract_draws(fit_main)

adj_draws <- adjust_params(
  draws, design = ct_model$design, onsets_flag = TRUE) 

adj_draws <- adj_draws %>% update_predictor_labels()

adj_draws <- add_baseline_to_draws(
  adj_draws, "Omicron (BA.1)", onsets_flag = TRUE)

adj_draws <- add_baseline_to_draws(
  adj_draws, "4 exposures", onsets_flag = TRUE)

adj_draws <- add_baseline_to_draws(
  adj_draws, "Symptomatic", onsets_flag = TRUE)

adj_draws <- add_baseline_to_draws(
  adj_draws, "Age: 35-49", onsets_flag = TRUE)

adj_draws <- adj_draws %>% 
  transform_to_model(., onsets_flag = TRUE) %>% 
  add_regressor_categories()

dt_inc_period <- adj_draws[,
  ip_draw := rlnorm(1, meanlog = inc_mean, sdlog = inc_sd), 
  by = c(".draw", "predictor", "regressor_category")]

factor_order <- c(
  "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Symptomatic", 
  "3 exposures", "4 exposures", "5+ exposures", "Age: 20-34", 
  "Age: 35-49", "Age: 50+")

dt_inc_period[, predictor := fct_relevel(predictor, factor_order)]

# Plotting incubation periods as ridge plots, using ggridges
p_ip_dist <- dt_inc_period[
  !is.na(regressor_category) & predictor != "Asymptomatic"] %>% 
  ggplot() + 
  geom_density_ridges(aes(x = ip_draw, y = predictor, fill = predictor)) + 
  # facet_wrap(~regressor_category) + 
  lims(x = c(0, 10)) + 
  theme_minimal() + 
  labs(x = "Time of symptom onset (days since exposure)", 
       y = "Probability",
       fill = "Predictor") +
  scale_y_discrete(limits = rev) + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

# Calculating effect sizes: medians and 95% credible intervals of overall 
# incubation periods
dt_inc_period[
  , `:=` (me = quantile(ip_draw, 0.5, na.rm = TRUE),
          lo = quantile(ip_draw, 0.025, na.rm = TRUE), 
          hi = quantile(ip_draw, 0.975, na.rm = TRUE)),
  by = c("predictor", "regressor_category")]

# Plotting the summarised incubation periods
p_ip_eff <- dt_inc_period[
  !is.na(regressor_category) & predictor != "Asymptomatic"] %>% 
  ggplot() + 
  geom_pointrange(aes(x = predictor,
                      y = me,
                      ymin = lo,
                      ymax = hi,
                      colour = predictor)) +
  # facet_wrap(~regressor_category) + 
  # lims(x = c(0, 10)) + 
  theme_minimal() + 
  labs(x = "Predictor", 
       y = "Time of symptom onset (days since exposure)",
       colour = "Predictor") + 
  scale_x_discrete(limits = rev) +
  coord_flip() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

# Putting the two panels together
p_ip <- grid.arrange(
  patchworkGrob(
    p_ip_dist + p_ip_eff + 
      plot_annotation(
        title = 'Incubation periods by covariates')),
  bottom = "Time of symptom onset (days since exposure)")

# Saving the data required to remake the figure
saveRDS(dt_inc_period, "outputs/plot_data/supplement/figure_S22.rds")

# Saving the final plot
ggsave("outputs/figures/pdfs/supplement/figure_S22.pdf",
       p_ip,
       width = 8,
       height = 8)


