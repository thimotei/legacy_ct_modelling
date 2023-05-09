# plotting the Ct value adjustment effect size posteriors

# script for making figure 3
library(ggplot2)
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
# library(facetscales)
library(gridExtra)
library(data.table)
library(ggExtra)

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup.R")

# either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
fit <- readRDS("outputs/fits/fit_full.rds")
draws <- extract_draws(fit)

# extracting beta coefficients from the fitted object
coeff_draws <- extract_coeffs(draws, 
                              design = adjustment_model$design)

# binding the set of draws from the fitted object corresponding to
# the baseline case
coeff_draws <- rbind(
  extract_coeffs(draws, 
                 design = adjustment_model$design)[,
                                                   predictor := "ORF1ab"
                 ],
  coeff_draws
)[, predictor := factor(predictor)]


coeff_draws_plot <- coeff_draws[variable %in% c("ct_shift", "ct_scale")]

# changing the names of the predictors and variables for the plot
coeff_draws_plot[predictor == "ct_typect_n_gene", predictor := "N gene"]
coeff_draws_plot[predictor == "ct_typect_s_gene", predictor := "S gene"]
coeff_draws_plot[predictor == "swab_typeVTM", predictor := "VTM swab"]
coeff_draws_plot[variable == "ct_scale", variable := "Ct value scaling adjustment"]
coeff_draws_plot[variable == "ct_shift", variable := "Ct value shift adjustment"]


# plotting the posteriors as a ridge plot
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

coeff_draws_plot_sum <- coeff_draws_plot[, 
  .(me = quantile(value, 0.5),
    lo = quantile(value, 0.025),
    hi = quantile(value, 0.975)), 
  by = c("variable", "predictor")]

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

p_adjustment <- plot_grid(p_adjustment_coeffs, p_adjustment_coeffs_sum, rel_widths = c(1, 0.7))

# saving the plot as a PDF
ggsave("outputs/figures/pdfs/ct_adjustment_posteriors.pdf",
       p_adjustment,
       width = 6,
       height = 6)

ggsave("outputs/figures/pngs/ct_adjustment_posteriors.png",
       p_adjustment,
       width = 6,
       height = 6)


