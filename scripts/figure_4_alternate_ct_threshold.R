# script for making figure 3
library(ggplot2)
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
library(facetscales)
library(data.table)
library(ggExtra)

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup.R")

# either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
fit <- readRDS("outputs/fits/fit_full.rds")
draws <- extract_draws(fit)

# only run is the script for figure 3 hasn't already been run,
# exactly the same command. Takes a long time — commented out by default
dt_pop_ct_draws <- extract_pop_ct_trajectories(fit,
                                               no_draws = 10000,
                                               onsets = TRUE)

dt_figure_4 <- figure_4_data(dt_pop_ct_draws, 
                             ct_threshold = 25)

# bivariate plot of incubation period and peak Ct value
figure_41 <- figure_4_panel(dt_figure_4,
                            regressor_category_arg = "VOC",
                            title_arg = "VOC",
                            factor_arg = c("Delta", 
                                           "Omicron (BA.1)",
                                           "Omicron (BA.2)")) +
  scale_colour_nejm()

figure_41_w_marginals <- add_marginals(figure_41)

figure_42 <- figure_4_panel(dt_figure_4,
                            regressor_category_arg = "Symptom status",
                            title_arg = "Symptom status",
                            factor_arg = "Symptomatic") +
  scale_colour_npg()

figure_42_w_marginals <- add_marginals(figure_42)

figure_43 <- figure_4_panel(dt_figure_4,
                            regressor_category_arg = "Number of exposures",
                            title_arg = "Number of exposures",
                            factor_arg = c("3 exposures",
                                           "4 exposures",
                                           "5+ exposures")) + 
  scale_colour_brewer(palette = "Set2")

figure_43_w_marginals <- add_marginals(figure_43)

figure_44 <- figure_4_panel(dt_figure_4,
                            regressor_category_arg = "Age",
                            title_arg = "Age",
                            factor_arg = c("Age: 20-34",
                                           "Age: 35-49",
                                           "Age: 50+")) + 
  scale_colour_aaas()

figure_44_w_marginals <- add_marginals(figure_44)

figure_4 <- plot_grid(figure_41_w_marginals,
                      figure_42_w_marginals,
                      figure_43_w_marginals,
                      figure_44_w_marginals,
                      ncol = 2)

ggsave("outputs/figures/pngs/figure_4_alt_ct_thres.png",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_4_alt_ct_thres.pdf",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")

