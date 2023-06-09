# script for making figure 3
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
library(facetscales)
library(data.table)

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup.R")

# Fit the model if it hasn't been run and no fit object exists locally
# fit_inf_lkj_prior <- epict(
#   obs,
#   ct_model = ct_model,
#   adjustment_model  = adjustment_model,
#   likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
#   onsets = TRUE, # include the symptom onset data and likelihood components
#   switch = FALSE, # include the longer wane part of the model
#   individual_variation = 0.025, # control the amount of individual-level variation
#   individual_correlation = 50, # no correlation between individual-level parameters
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
# fit_inf_lkj_prior$save_object("outputs/fits/fit_inf_lkj_prior.rds")

# reading in fit object, if already run
fit_inf_lkj_prior <- readRDS("outputs/fits/fit_inf_lkj_prior.rds")

# either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
draws <- extract_draws(fit_inf_lkj_prior)

# the next command takes a long time. A high number of draws was used 
# (no_draws) for the figures but isn't necessary for quicker plots
pop_ct_draws_trim <- extract_pop_ct_trajectories(fit_inf_lkj_prior,
                                                 no_draws = 10000,
                                                 onsets = TRUE)
# summarising Ct trajectories
pop_ct_draws_sum <- summarise_ct_traj(pop_ct_draws_trim)

# summarising effect sizes and converting to natural units
effect_size_summary_natural <- summarise_effect_sizes_natural(draws)

# VOC panel
p3_11_inf_lkj_prior <- plot_ct_trajectory_panel(
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

p3_12_inf_lkj_prior <- plot_effect_panel(
  effect_size_summary_natural, 
  regressor_category_arg = "VOC", 
  baseline_arg = "Omicron (BA.1)",
  factor_arg = c("Delta",
                 "Omicron (BA.1)",
                 "Omicron (BA.2)")) +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Symptoms panel
p3_21_inf_lkj_prior <- plot_ct_trajectory_panel(
  pop_ct_draws_trim,
  pop_ct_draws_sum,
  regressor_category_arg = "Symptom status",
  no_draws = 500,
  factor_arg = "Symptomatic",
  title_arg = "Symptom status",
  alpha_arg = 0.15) +
  scale_colour_npg() + 
  scale_fill_npg()

p3_22_inf_lkj_prior <- plot_effect_panel(
  effect_size_summary_natural, 
  regressor_category_arg = "Symptom status",
  baseline_arg = "Symptomatic",
  c("Symptomatic", "Asymptomatic")) +
  scale_colour_npg() 

# Number of exposures panel
p3_31_inf_lkj_prior <- plot_ct_trajectory_panel(
  pop_ct_draws_trim,
  pop_ct_draws_sum,
  regressor_category_arg = "Number of exposures",
  no_draws = 500,
  factor_arg = "3 exposures",
  title_arg = "Number of exposures",
  alpha_arg = 0.15) +
  scale_colour_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") 

p3_32_inf_lkj_prior <- plot_effect_panel(
  effect_size_summary_natural, 
  regressor_category_arg = "Number of exposures",
  baseline_arg = "4 exposures",
  c("3 exposures", "4 exposures",
    "5+ exposures")) +
  scale_fill_brewer(palette = "Set2") 

# Age panel
p3_41_inf_lkj_prior <- plot_ct_trajectory_panel(
  pop_ct_draws_trim,
  pop_ct_draws_sum,
  regressor_category_arg = "Age",
  no_draws = 500,
  factor_arg = "Age: 20-34",
  title_arg = "Age",
  alpha_arg = 0.15) +
  scale_colour_aaas() +
  scale_fill_aaas()

p3_42_inf_lkj_prior <- plot_effect_panel(
  effect_size_summary_natural,
  regressor_category_arg = "Age",
  baseline_arg = "Age: 35-49",
  c("Age: 20-34", "Age: 35-49",
    "Age: 50+")) +
  scale_colour_aaas()

# plotting the trajectory panels and effect size panels together
p3_1_inf_lkj_prior <- figure_3_subpanel(
  p3_11_inf_lkj_prior, 
  p3_12_inf_lkj_prior, 
  rel_widths_arg = c(2, 1))

p3_2_inf_lkj_prior <- figure_3_subpanel(
  p3_21_inf_lkj_prior, 
  p3_22_inf_lkj_prior,
  rel_widths_arg = c(2, 1))

p3_3_inf_lkj_prior <- figure_3_subpanel(
  p3_31_inf_lkj_prior, 
  p3_32_inf_lkj_prior,
  rel_widths_arg = c(2, 1))

p3_4_inf_lkj_prior <- figure_3_subpanel(
  p3_41_inf_lkj_prior, 
  p3_42_inf_lkj_prior,
  rel_widths_arg = c(2, 1))

# plotting all 4 regressor category panels together
p31_inf_lkj_prior <- plot_grid(
  p3_1_inf_lkj_prior,
  p3_2_inf_lkj_prior,
  p3_3_inf_lkj_prior,
  p3_4_inf_lkj_prior, nrow = 2)

# saving final figure
ggsave("outputs/figures/pngs/figure_3_inf_lkj_prior.png",
       p31_inf_lkj_prior,
       width = 12,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_3_inf_lkj_prior.pdf",
       p31_inf_lkj_prior,
       width = 12,
       height = 10,
       bg = "white")

