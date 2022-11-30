# script for making figure 3
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
library(facetscales)
library(data.table)

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup_no_covariates.R")

# Fit the model if it hasn't been run and no fit object exists locally
# fit_no_covariates <- epict(
#   obs,
#   ct_model = ct_model,
#   adjustment_model  = adjustment_model,
#   likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
#   onsets = TRUE, # include the symptom onset data and likelihood components
#   switch = FALSE, # include the longer wane part of the model
#   individual_variation = 0.025, # control the amount of
#   individual_correlation = 2,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 2000,
#   adapt_delta = 0.95,
#   max_treedepth = 12,
#   output_loglik = FALSE
# )

# saves fitted object - usually results in a very large file (> 100Mbs)
# if the full model is fit
# fit_no_covariates$save_object("outputs/fits/fit_no_covariates.rds")

# reading in fit object, if already run
fit_no_covariates <- readRDS("outputs/fits/fit_no_covariates.rds")

# either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
draws <- extract_draws(fit_no_covariates)

# the next command takes a long time. A high number of draws was used 
# (no_draws) for the figures but isn't necessary for quicker plots
pop_ct_draws_trim <- extract_pop_ct_trajectories(
  fit_no_covariates,
  no_draws = 10000,
  separate_baseline_covariates = FALSE,
  other_covariates = FALSE,
  onsets = TRUE)

# summarising Ct trajectories
pop_ct_draws_sum <- summarise_ct_traj(pop_ct_draws_trim)

# summarising effect sizes and converting to natural units
effect_size_summary_natural <- summarise_effect_sizes_natural(
  draws, add_baseline_flag = FALSE)

# VOC panel
p3_11_no_covariates <- plot_ct_trajectory_panel(
  pop_ct_draws_cens = pop_ct_draws_trim,
  pop_ct_draws_summ = pop_ct_draws_sum,
  regressor_category_arg = "VOC", 
  no_draws = 500,
  factor_arg = c("Delta",
                 "Omicron (BA.1)",
                 "Omicron (BA.2)"),
  title_arg = "VOC",
  alpha_arg = 0.15) +
  scale_colour_nejm() + 
  scale_fill_nejm()

p3_12_no_covariates <- plot_effect_panel(
  effect_size_summary_natural, 
  regressor_category_arg = "VOC", 
  baseline_arg = "Omicron (BA.1)",
  factor_arg = c("Delta",
                 "Omicron (BA.1)",
                 "Omicron (BA.2)")) +
  scale_colour_nejm() + 
  scale_fill_nejm()

# plotting the trajectory panels and effect size panels together
p3_1_no_covariates <- figure_3_subpanel(
  p3_11_no_covariates, 
  p3_12_no_covariates, 
  rel_widths_arg = c(2, 1))

# saving final figure
ggsave("outputs/figures/pngs/figure_3_no_covariates.png",
       p3_1_no_covariates,
       width = 8,
       height = 6,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_3_no_covariates.pdf",
       p3_1_no_covariates,
       width = 8,
       height = 6,
       bg = "white")

