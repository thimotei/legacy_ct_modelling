#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(stringr)
library(purrr)
library(lubridate)
library(patchwork)
library(tidyr)
library(here)

# loading all functions in package directory
files <- list.files("R", "*.R", full.names = TRUE)
walk(files, source)

# load in data processed in scripts/process-data.R
dt_clean <- fread(here("data/processed-data.csv"))

# Do additional processing to filter for the desired number of swabs
dt_2_tests <- subset_data(
  dt_clean, voc = c("Delta", "Omicron"), no_pos_swabs = 2
)

# Plot the raw data
p1_raw <- plot_obs_ct(dt_2_tests)

# Specify the CT model design matrix
ct_model <- subject_design(
  ~ 1 + VOC,
  data = dt_2_tests,
  params = c("t_p", "t_s", "c_p", "c_s", "inc_mean"),
  preds_sd = 0.1
)

# Translate data and model specification to stan format
stan_data <- data_to_stan(
  dt_2_tests, ct_model = ct_model,
  likelihood = TRUE, onsets = TRUE
)

# Compile model
mod <- cmdstan_model("stan/ct_trajectory_model.stan", include_paths = "stan")

# Fit
fit <- mod$sample(
  data = stan_data,
  init = stan_inits(stan_data),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 15
)

# Print parameter summaries
summarise_pop_pp(fit)

summarise_coeff_pp(fit)

# Extract and plot posterior predictions
ct_draws <- extract_ct_trajectories(fit)

ct_summary <- summarise_draws(
  copy(ct_draws)[,
    time_since_first_pos := as.integer(time_since_first_pos)
    ],
  by = c("id", "time_since_first_pos")
)

ct_pp <- extract_posterior_predictions(fit, dt_2_tests)
ct_pp <- summarise_draws(
  ct_pp[, value := sim_ct], by = c("id", "t", "pcr_res", "obs")
)

# Plotting summaries of fitted trajectories against simulated data
pp_plot <- plot_obs_ct(
  dt_2_tests, ct_draws[iteration <= 10], ct_pp, traj_alpha = 0.05
)

ggsave("outputs/figures/pp.png", pp_plot, height = 16, width = 16)

# Extract and plot population level posterior predictions for the CT model
draws <- extract_draws(fit)

pop_pp <- plot_ct_summary(
  draws, time_range = seq(0, 60, by = 0.01), samples = 100, by = c()
)

ggsave(
  "outputs/figures/population_ct_pp.png",
  pop_pp, width = 8, height = 8,
)
