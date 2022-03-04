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
devtools::load_all()

# load in data processed in scripts/process-data.R
dt_clean <- readRDS(here("data/processed-data.rds"))

# Do additional processing to filter for the desired number of swabs
# per positive episode
dt_2_tests <- subset_data(dt_clean, no_pos_swabs = 2)

# Plot the raw data
p1_raw <- plot_obs_ct(dt_2_tests) +
   facet_wrap(vars(id, VOC))

# Specify which params adjusting for (see params_avail_to_adjust() for options)
# Here all available options (can also specify this using "all")
adj_params <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")

# Specify the CT model design matrix
ct_model <- subject_design(
  ~ 1 + VOC + symptoms + no_vaccines,
  data = dt_2_tests,
  params = adj_params,
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
summarise_coeff_pp(fit, params = adj_params, exponentiate = TRUE)

# Extract and plot posterior predictions
pp_plot <- plot_pp_from_fit(
  fit, obs = dt_2_tests, samples = 50, alpha = 0.025
) + 
  facet_wrap(vars(id))

ggsave("outputs/figures/pp.png", pp_plot, height = 16, width = 16)

# Extract and plot population level posterior predictions for the CT model
draws <- extract_draws(fit)

parameter_pp <- plot_summary(draws)

ggsave(
  "outputs/figures/parameter_pp.png",
  parameter_pp, width = 12, height = 8,
)

# Extract effect sizes and make a summary plot
eff_plot <- draws %>%
  summarise_effects(design = ct_model$design) %>%
  update_variable_labels(reverse = TRUE) %>%
  plot_effects() +
  facet_wrap(vars(preds))

ggsave(
  "outputs/figures/effects_summary.png",
  eff_plot, width = 16, height = 16,
)