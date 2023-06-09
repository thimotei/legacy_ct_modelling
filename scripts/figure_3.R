# script for making figure 3
library(cowplot)
library(ggridges)
library(lemon)
library(scico)
library(ggsci)
library(data.table)

# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup.R")

# either run inference script, which saves a fit object, or load a saved fit
# object from a previous inference
fit <- readRDS("outputs/fits/fit_full.rds")
draws <- extract_draws(fit)

# the next command takes a long time. A high number of draws was used 
# (no_draws) for the figures
dt_pop_ct_draws <- extract_pop_ct_trajectories(fit,
                                               no_draws = 10000,
                                               onsets = TRUE)

# summarising Ct trajectories
pop_ct_draws_sum <- summarise_ct_traj(dt_pop_ct_draws, pop_flag = TRUE)

# summarising effect sizes and converting to natural units
effect_size_summary_natural <- summarise_effect_sizes_natural(draws)

# VOC panel
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

# Symptoms panel
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

# Number of exposures panel
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

# Age panel
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

p3_1 <- figure_3_subpanel(p3_11, p3_12, rel_widths_arg = c(2, 1))
p3_2 <- figure_3_subpanel(p3_21, p3_22, rel_widths_arg = c(2, 1))
p3_3 <- figure_3_subpanel(p3_31, p3_32, rel_widths_arg = c(2, 1))
p3_4 <- figure_3_subpanel(p3_41, p3_42, rel_widths_arg = c(2, 1))

p31 <- plot_grid(p3_1, p3_2, p3_3, p3_4, nrow = 2)

ggsave("outputs/figures/pngs/figure_3_new.png",
       p31,
       width = 12,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_3_new.pdf",
       p31,
       width = 12,
       height = 10,
       bg = "white")

