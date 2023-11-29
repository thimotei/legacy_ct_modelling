#--- Figure S16
#--- Plotting the Ct value adjustment effect size posteriors

# Sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/voc.R")

# Run the script found at /scripts/inference/loo_comparison/no_covariates.R 
# if you do not have a saved .rds object containing the cmdstanr fit object
# after fitting the model to the dataset

# Reading in fit object, if already run
fit_voc <- readRDS("outputs/fits/voc.rds")

# Either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
draws <- extract_draws(fit_voc)

# The next command takes a long time. A high number of draws was used 
# (no_draws) for the figures but isn't necessary for quicker plots
dt_pop_ct_draws <- extract_pop_ct_trajectories(
  fit_voc, no_draws = 2000, 
  separate_baseline_covariates = FALSE,
  other_covariates = FALSE,
  onsets = TRUE)

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
p3_11_no_covariates <- plot_ct_trajectory_panel(
  pop_ct_draws_cens = dt_pop_ct_draws,
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

# Saving all three required data.tables to build the overall plot
saveRDS(dt_pop_ct_draws, "outputs/plot_data/supplement/figure_S16_traj.rds")
saveRDS(pop_ct_draws_sum, "outputs/plot_data/supplement/figure_S16_sum.rds")
saveRDS(effect_size_summary_natural, "outputs/plot_data/supplement/figure_S16_eff.rds")

# Saving final figure
ggsave("outputs/figures/pngs/supplement/figure_S16.png",
       p3_1_no_covariates,
       width = 12,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/figure_S16.pdf",
       p3_1_no_covariates,
       width = 12,
       height = 10,
       bg = "white")


