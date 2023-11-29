#--- Figure S23
#--- Plotting the Ct value adjustment effect size posteriors

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# Run the script found at /scripts/inference/main.R if you do not have a saved
# .rds object containing the cmdstanr fit object after fitting the model 
# to the dataset

# However, if you have the fit object already, just load in the data
fit_main <- readRDS("outputs/fits/main.rds")

# Extracting the posterior draws from the fit object
draws <- extract_draws(fit_main)

# Extracting just the beta coefficients from the fit object
coeff_draws <- extract_coeffs(
  draws, design = adjustment_model$design)

# Binding the set of draws from the fitted object corresponding to
# the baseline case
coeff_draws <- rbind(
  extract_coeffs(
    draws, design = adjustment_model$design)[, predictor := "ORF1ab"], 
  coeff_draws
)[, predictor := factor(predictor)]

# Filtering for the Ct adjustment parmaeters
coeff_draws_plot <- coeff_draws[variable %in% c("ct_shift", "ct_scale")]

# Changing the names of the predictors and variables for the plot
coeff_draws_plot[predictor == "ct_typect_n_gene", predictor := "N gene"]
coeff_draws_plot[predictor == "ct_typect_s_gene", predictor := "S gene"]
coeff_draws_plot[predictor == "swab_typeVTM", predictor := "VTM swab"]
coeff_draws_plot[variable == "ct_scale", variable := "Ct value scaling adjustment"]
coeff_draws_plot[variable == "ct_shift", variable := "Ct value shift adjustment"]

# Plotting the posteriors as a ridge plot
p_adjustment_coeffs <- coeff_draws_plot |> 
  ggplot() +
  geom_density_ridges(aes(x = value,
                          y = predictor, 
                          fill = predictor),
                      alpha = 0.5) + 
  facet_wrap(~variable, scales = "free", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "value",
       y = "Gene target/swab type",
       fill = "Adjusted parameter") + 
  scale_y_discrete(limits = rev)

# Summarising the posteriors by calculating the median
# and 95% credible interval
coeff_draws_plot_sum <- coeff_draws_plot[, 
  .(me = quantile(value, 0.5),
    lo = quantile(value, 0.025),
    hi = quantile(value, 0.975)), 
  by = c("variable", "predictor")]

# Plotting the median and credible interval as a pointrange
p_adjustment_coeffs_sum <- coeff_draws_plot_sum |> 
  ggplot() +
  geom_pointrange(aes(x = predictor,
                      y = me,
                      ymin = lo,
                      ymax = hi,
                      colour = predictor)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  # geom_blank(data = blank_data_plot,
  #            aes(x = predictor, y = y)) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  facet_rep_wrap(~variable,
                 scales = "free",
                 ncol = 1,
                 strip.position = "bottom") + 
  # scale_color_aaas() +
  custom_plot_theme() +
  labs(y = "Density") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# Plotting the ridgeplot next to the credible intervals
p_adjustment <- plot_grid(
  p_adjustment_coeffs, 
  p_adjustment_coeffs_sum, 
  rel_widths = c(1, 0.7))

saveRDS(coeff_draws_plot, "outputs/plot_data/supplement/figure_S23.rds")
saveRDS(coeff_draws_plot_sum, "outputs/plot_data/supplement/figure_S23_summary.rds")

# Saving the final plot
ggsave("outputs/figures/pdfs/supplement/figure_S23.pdf",
       p_adjustment,
       width = 6,
       height = 6,
       bg = "white")

ggsave("outputs/figures/pngs/supplement/figure_S23.png",
       p_adjustment,
       width = 6,
       height = 6,
       bg = "white")


