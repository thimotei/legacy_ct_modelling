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
# draws <- extract_draws(fit)

# only run is the script for figure 3 hasn't already been run,
# exactly the same command. Takes a long time — commented out by default
dt_pop_ct_draws <- extract_pop_ct_trajectories(fit,
                                               no_draws = 10000,
                                               onsets = FALSE)

dt_figure_4 <- munge_figure_4_data(draws,
                                   dt_pop_ct_draws,
                                   ct_threshold_arg = 20)

# bivariate plot of incubation period and peak Ct value
figure_41_tmp <- figure_4_panel(dt_figure_4,
                                regressor_category_arg = "VOC",
                                title_arg = "VOC",
                                factor_arg = c("Delta", 
                                               "Omicron (BA.1)",
                                               "Omicron (BA.2)")) +
  scale_colour_nejm()

figure_41 <- add_marginals(figure_41_tmp)

figure_42_tmp <- figure_4_panel(dt_figure_4,
                                regressor_category_arg = "Symptom status",
                                title_arg = "Symptom status",
                                factor_arg = "Symptomatic") +
  scale_colour_npg()

figure_42 <- add_marginals(figure_42_tmp)

figure_43_tmp <- figure_4_panel(dt_figure_4,
                                regressor_category_arg = "Number of exposures",
                                title_arg = "Number of exposures",
                                factor_arg = c("3 exposures",
                                               "4 exposures",
                                               "5+ exposures")) + 
  scale_colour_brewer(palette = "Set2")

figure_43 <- add_marginals(figure_43_tmp)

figure_44_tmp <- figure_4_panel(dt_figure_4,
                                regressor_category_arg = "Age",
                                title_arg = "Age",
                                factor_arg = c("Age: 20-34",
                                               "Age: 35-49",
                                               "Age: 50+")) + 
  scale_colour_aaas()

figure_44 <- add_marginals(figure_44_tmp)

figure_4 <- plot_grid(figure_41, figure_42, figure_43, figure_44)

ggsave("outputs/figures/pngs/figure_4_new.png",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_4_new.pdf",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")
