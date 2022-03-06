#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(cowplot)
library(stringr)
library(purrr)

# loading all functions in package directory
devtools::load_all() 

# Set up simulations for 20 individuals
# To be lazy here we are assuming that onsets are not available
# To update need another modelling step for onsets
obs <- list(
  P = 20,
  any_onsets = 1,
  onset_time = rep(0, 20),
  c_lod = 40,
  lmean = get_inc_period()$inc_mean_p,
  lsd = get_inc_period()$inc_sd_p,
  swab_types = 0,
  preds = 0
)

# Simulate from the centre of the prior for all parameters
# based on initial conditions used for the stan model.
ct_sample <- simulate_obs(
  obs = obs,
  parameters = stan_inits(obs)(),
  time_range = -1:30,
  sample_density = 4:12
)

# plot of subset of data
plot_obs_ct(ct_sample) +
  facet_wrap(vars(id)) +

# compiling model
mod <- cmdstan_model(
  "stan/ct_trajectory_model.stan",
  include_paths = "stan",
  stanc_options = list("O1")
)

sim_stan_data <- data_to_stan(ct_sample, onset = TRUE)

# fitting the model - not very quick, as many iterations hit the
# max_tree_depth at the moment
fit_sim <- mod$sample(
  data = sim_stan_data,
  init = stan_inits(sim_stan_data),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Extract and plot posterior predictions
sim_pp_plot <- plot_pp_from_fit(
  fit_sim, obs = ct_sample, samples = 50, alpha = 0.025
) +
  facet_wrap(vars(id)) +

ggsave("outputs/figures/sim_pp.png", sim_pp_plot, height = 10, width = 10)

# Extract and plot population level posterior predictions for the CT model
sim_draws <- extract_draws(fit_sim)

sim_parameter_pp <- plot_summary(sim_draws)

ggsave(
  "outputs/figures/sim_parameter_pp.png",
  sim_parameter_pp, width = 12, height = 8,
)


ggsave(
  "outputs/figures/sim_population_ct_pp.png",
  sim_pop_pp, width = 8, height = 8,
)
