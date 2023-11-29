#--- Script for making figure 3
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
library(data.table)

# Sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# Run the script found at /scripts/inference/main.R if you do not have a saved
# .rds object containing the cmdstanr fit object after fitting the model 
# to the dataset

# However, if you have the fit object already, just load in the data
fit_main <- readRDS("outputs/fits/main.rds")

# Extracting the posterior draws from the fit object
draws <- extract_draws(fit_main)

#--- Computing Ct trajectories from the population-level posteriors

# This function takes the draws, puts them back to a natural scale
# and simulates the underlying model for each set of population-level
# process parameters. It takes a while to run
dt_pop_ct_draws <- extract_pop_ct_trajectories(
  fit_main, no_draws = 2000, onsets = TRUE)

# Summarising the many Ct trajectories, computing the median and 
# 95% credible interval for each covariate
pop_ct_draws_sum <- summarise_ct_traj(dt_pop_ct_draws, pop_flag = TRUE)

# Calculating time trajectories hit the lab-based LOD (not the time our latent
# LOD, which is what the parameter that comes directly out of the inference
# represents) 
dt_t_lod_sum <- figure_3_t_lod_sum(dt_pop_ct_draws, ct_threshold = 37)
dt_t_lod_raw <- figure_3_t_lod(dt_pop_ct_draws, ct_threshold = 37)

# Adding the baseline category to the new t_lod data.table
dt_t_lod_baseline <- add_baseline_effect_sizes(dt_t_lod_raw)

# Summarising effect sizes and converting to natural units
effect_size_summary_natural_tmp <- summarise_effect_sizes_natural(
  draws, 
  ct_model = ct_model,
  add_baseline_flag = TRUE,
  onsets_flag = FALSE)

# Removing the time of the latent LOD from the effect size panel, as 
# we are interested in the lab-based LOD
effect_size_summary_natural_tmp <- effect_size_summary_natural_tmp[
  variable != "Time until PCR-"]

# Combining the lab and latent time of LODs
effect_size_summary_natural <- rbind(
  effect_size_summary_natural_tmp, dt_t_lod_sum,
  dt_t_lod_baseline)

#--- Plotting the trajectories and effect sizes together, one panel at a time
#--- VOC panel
p3_11 <- plot_ct_trajectory_panel(pop_ct_draws_cens = dt_pop_ct_draws,
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

p3_12 <- plot_effect_panel(effect_size_summary_natural, 
                           regressor_category_arg = "VOC", 
                           baseline_arg = "Omicron (BA.1)",
                           factor_arg = c("Delta",
                                          "Omicron (BA.1)",
                                          "Omicron (BA.2)")) +
  scale_colour_nejm() + 
  scale_fill_nejm()

#--- Symptoms panel
p3_21 <- plot_ct_trajectory_panel(dt_pop_ct_draws,
                                  pop_ct_draws_sum,
                                  regressor_category_arg = "Symptom status",
                                  no_draws = 500,
                                  factor_arg = "Symptomatic",
                                  title_arg = "Symptom status",
                                  alpha_arg = 0.15) +
  scale_colour_npg() + 
  scale_fill_npg()

p3_22 <- plot_effect_panel(effect_size_summary_natural, 
                           regressor_category_arg = "Symptom status",
                           baseline_arg = "Symptomatic",
                           c("Symptomatic", "Asymptomatic")) +
  scale_colour_npg() 

#--- Number of exposures panel
p3_31 <- plot_ct_trajectory_panel(dt_pop_ct_draws,
                                  pop_ct_draws_sum,
                                  regressor_category_arg = "Number of exposures",
                                  no_draws = 500,
                                  factor_arg = "3 exposures",
                                  title_arg = "Number of exposures",
                                  alpha_arg = 0.15) +
  scale_colour_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") 

p3_32 <- plot_effect_panel(effect_size_summary_natural, 
                           regressor_category_arg = "Number of exposures",
                           baseline_arg = "4 exposures",
                           c("3 exposures", "4 exposures",
                             "5+ exposures")) +
  scale_fill_brewer(palette = "Set2") 

#--- Age panel
p3_41 <- plot_ct_trajectory_panel(dt_pop_ct_draws,
                                  pop_ct_draws_sum,
                                  regressor_category_arg = "Age",
                                  no_draws = 500,
                                  factor_arg = "Age: 20-34",
                                  title_arg = "Age",
                                  alpha_arg = 0.15) +
  scale_colour_aaas() +
  scale_fill_aaas()

p3_42 <- plot_effect_panel(effect_size_summary_natural,
                           regressor_category_arg = "Age",
                           baseline_arg = "Age: 35-49",
                           c("Age: 20-34", "Age: 35-49",
                             "Age: 50+")) +
  scale_colour_aaas()

# Putting the trajectory subpanels together with the summary subpanels
p3_1 <- figure_3_subpanel(p3_11, p3_12, rel_widths_arg = c(2, 1))
p3_2 <- figure_3_subpanel(p3_21, p3_22, rel_widths_arg = c(2, 1))
p3_3 <- figure_3_subpanel(p3_31, p3_32, rel_widths_arg = c(2, 1))
p3_4 <- figure_3_subpanel(p3_41, p3_42, rel_widths_arg = c(2, 1))

# Putting all four panels together
p31 <- plot_grid(p3_1, p3_2, p3_3, p3_4, nrow = 2)

# Saving the final figure
ggsave("outputs/figures/pngs/figure_3.png",
       p31,
       width = 12,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_3.pdf",
       p31,
       width = 12,
       height = 10,
       bg = "white")
