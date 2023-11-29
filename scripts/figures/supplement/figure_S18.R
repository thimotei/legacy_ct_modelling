#--- Figure S18
#--- Plotting the Ct value adjustment effect size posteriors

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# Run the script found at /scripts/inference/loo_comparison/inf_lkj_prior.R if 
# you do not have a saved rds object containing the cmdstanr fit object after 
# fitting the model to the dataset

# However, if you have the fit object already, just load in the data
# reading in fit object, if already run
fit_inf_lkj_prior <- readRDS("outputs/fits/inf_lkj_prior.rds")

# Extract posterior samples from fit object
draws <- extract_draws(fit_inf_lkj_prior)

# The next command takes a long time. A high number of draws was used 
# (no_draws) for the figures but isn't necessary for quicker plots
dt_pop_ct_draws <- extract_pop_ct_trajectories(
  fit_inf_lkj_prior, no_draws = 2000, onsets = TRUE)

# Summarising Ct trajectories
pop_ct_draws_sum <- summarise_ct_traj(dt_pop_ct_draws)

#--- summarising the population-level posteriors and representing them
#--- as effect sizes relative to the baseline set of individuals 

# calculating time trajectories hit the lab-based LOD (not the time our latent
# LOD, which is what the parameter that comes directly out of the inference
# represents) 
dt_t_lod_sum <- figure_3_t_lod_sum(dt_pop_ct_draws, ct_threshold = 37)
dt_t_lod_raw <- figure_3_t_lod(dt_pop_ct_draws, ct_threshold = 37)

# adding the baseline category to the new t_lod data.table
dt_t_lod_baseline <- add_baseline_effect_sizes(dt_t_lod_raw)

# summarising effect sizes and converting to natural units
effect_size_summary_natural_tmp <- summarise_effect_sizes_natural(
  draws, 
  ct_model = ct_model,
  add_baseline_flag = TRUE,
  onsets_flag = FALSE)

# removing the time of the latent LOD from the effect size panel, as 
# we are interested in the lab-based LOD
effect_size_summary_natural_tmp <- effect_size_summary_natural_tmp[
  variable != "Time until PCR-"]

# combining the lab and latent time of LODs
effect_size_summary_natural <- rbind(
  effect_size_summary_natural_tmp, dt_t_lod_sum,
  dt_t_lod_baseline)

#--- Plotting fits
# VOC panel
p3_11_inf_lkj_prior <- plot_ct_trajectory_panel(
  pop_ct_draws_cens = dt_pop_ct_draws,
  pop_ct_draws_summ = pop_ct_draws_sum,
  regressor_category_arg = "VOC", 
  no_draws = 500,
  factor_arg = c(
    "Delta", "Omicron (BA.1)", "Omicron (BA.2)"),
  title_arg = "VOC",
  alpha_arg = 0.15) +
  scale_colour_nejm() + 
  scale_fill_nejm()

p3_12_inf_lkj_prior <- plot_effect_panel(
  effect_size_summary_natural, 
  regressor_category_arg = "VOC", 
  baseline_arg = "Omicron (BA.1)",
  factor_arg = c(
    "Delta", "Omicron (BA.1)","Omicron (BA.2)")) +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Symptoms panel
p3_21_inf_lkj_prior <- plot_ct_trajectory_panel(
  dt_pop_ct_draws,
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
  dt_pop_ct_draws,
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
  c("3 exposures", "4 exposures", "5+ exposures")) +
  scale_fill_brewer(palette = "Set2") 

# Age panel
p3_41_inf_lkj_prior <- plot_ct_trajectory_panel(
  dt_pop_ct_draws,
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
  c("Age: 20-34", "Age: 35-49", "Age: 50+")) +
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

# Saving all three required data.tables to build the overall plot
saveRDS(dt_pop_ct_draws, "outputs/plot_data/supplement/figure_S18_traj.rds")
saveRDS(pop_ct_draws_sum, "outputs/plot_data/supplement/figure_S18_sum.rds")
saveRDS(effect_size_summary_natural, "outputs/plot_data/supplement/figure_S18_eff.rds")

# Saving final figure
ggsave("outputs/figures/pngs/supplement/figure_S18.png",
       p31_inf_lkj_prior,
       width = 12,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/figure_S18.pdf",
       p31_inf_lkj_prior,
       width = 12,
       height = 10,
       bg = "white")


