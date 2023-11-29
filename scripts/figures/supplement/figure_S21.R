library(ggExtra)

# Figure S21

# Sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# Run the script found at /scripts/inference/main.R if you do not have a saved
# .rds object containing the cmdstanr fit object after fitting the model 
# to the dataset

# However, if you have the fit object already, just load in the data
fit_main <- readRDS("outputs/fits/main.rds")

# Extracting the posterior draws from the fit object
draws <- extract_draws(fit_main)

# Only run is the script for figure 3 hasn't already been run. It has
# the same command to simulate many trajectories. It takes a long time
# to run, which is why it is commented out by default
# dt_pop_ct_draws <- extract_pop_ct_trajectories(
#   fit_main, no_draws = 10000, onsets = TRUE)

# Using the Ct trajectories to calculate the incubation period draws and
# peak Ct values
dt_figure_4 <- figure_4_data(dt_pop_ct_draws, ct_threshold = 25)

#--- Bivariate plots of incubation period samples and peak Ct value samples,
#--- one panel at a time for each covariate category, much the same as 
#--- the Figure 2 script

# Bivariate plot of incubation period and peak Ct value
# VOC panel
figure_41 <- figure_4_panel(
  dt_figure_4, regressor_category_arg = "VOC", 
  title_arg = "VOC", factor_arg = c(
    "Delta", "Omicron (BA.1)", "Omicron (BA.2)")) +
  scale_colour_nejm()

# Adding marginal distributions - using ggExtra package
figure_41_w_marginals <- add_marginals(figure_41)

# Symptom status panel
figure_42 <- figure_4_panel(
  dt_figure_4, regressor_category_arg = "Symptom status",
  title_arg = "Symptom status", factor_arg = "Symptomatic") +
  scale_colour_npg()

# Adding marginal distributions - using ggExtra package
figure_42_w_marginals <- add_marginals(figure_42)

# Number of exposures panel
figure_43 <- figure_4_panel(
  dt_figure_4, regressor_category_arg = "Number of exposures",
  title_arg = "Number of exposures", factor_arg = c(
    "3 exposures", "4 exposures", "5+ exposures")) + 
  scale_colour_brewer(palette = "Set2")

# Adding marginal distributions - using ggExtra package
figure_43_w_marginals <- add_marginals(figure_43)

# Adding marginal distributions - using ggExtra package
figure_44 <- figure_4_panel(
  dt_figure_4, regressor_category_arg = "Age", 
  title_arg = "Age", factor_arg = c(
    "Age: 20-34", "Age: 35-49", "Age: 50+")) + 
  scale_colour_aaas()

# Adding marginal distributions - using ggExtra package
figure_44_w_marginals <- add_marginals(figure_44)

# Putting the four panels together
figure_4 <- plot_grid(
  figure_41_w_marginals, figure_42_w_marginals, 
  figure_43_w_marginals, figure_44_w_marginals,
  ncol = 2)

# Saving all three required data.tables to build the overall plot
saveRDS(dt_figure_4, "outputs/plot_data/supplement/figure_S21.rds")

ggsave("outputs/figures/pdfs/supplement/figure_S21.pdf",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")

ggsave("outputs/figures/pngs/supplement/figure_S21.png",
       figure_4,
       width = 11,
       height = 8,
       bg = "white")

