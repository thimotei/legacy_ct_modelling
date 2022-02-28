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

# Specify which params adjusting for (see params_avail_to_adjust() for options)
adj_params <- "all"

# Specify the CT model design matrix
ct_model <- subject_design(
  ~ 1 + VOC,
  data = dt_2_tests,
  params = "all",
  preds_sd = 0.1
)

# Translate data and model specification to stan format
stan_data <- data_to_stan(
  dt_2_tests,
  ct_model = ct_model,
  likelihood = TRUE,
  onsets = TRUE
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

# Population level parameter summary
summarise_pop_pp(fit)

# Population level adjustment summary
summarise_coeff_pp(fit, params = adj_params)

# Extract and plot posterior predictions
pp_plot <- plot_pp_from_fit(fit, obs = dt_2_tests, samples = 10, alpha = 0.05)

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
