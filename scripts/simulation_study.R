#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(cowplot)
library(stringr)
library(purrr)

# loading all functions in package directory
files <- list.files("R", "*.R", full.names = TRUE)
walk(files, source)

# Set up simulations for 20 individuals
# To be lazy here we are assuming that onsets are not available
# To update need another modelling step for onsets
obs <- list(
  P = 20,
  any_onsets = 1,
  onset_time = rep(0, 20),
  c_lod = 40,
  lmean = get_inc_period()$inc_mean_p,
  lsd = get_inc_period()$inc_sd_p
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
plot_obs_ct(ct_sample)

# compiling model
mod <- cmdstan_model("stan/ct_trajectory_model.stan", include_paths = "stan")

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

# extracting Ct fits. Bit slow as it is at the moment
ct_draws <- extract_ct_trajectories(fit_sim)

# summarising trajectories using median and 95% CrI
ct_summary <- summarise_draws(
  copy(ct_draws)[,
    time_since_first_pos := as.integer(time_since_first_pos)
    ],
  by = c("id", "time_since_first_pos")
)

# extract posterior CT predictons and  summarise
ct_pp <- extract_posterior_predictions(fit_sim, ct_sample)
ct_pp <- summarise_draws(ct_pp[, value := sim_ct], by = c("id", "t", "pcr_res"))

# plotting summaries of fitted trajectories against simulated data
sim_pp_plot <- plot_obs_ct(
  ct_sample, ct_draws[iteration <= 10], ct_pp, traj_alpha = 0.05
)
ggsave("outputs/figures/sim_pp.png", sim_pp_plot, height = 10, width = 10)