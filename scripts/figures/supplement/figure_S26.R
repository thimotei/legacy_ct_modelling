#--- Simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(cowplot)
library(stringr)
library(purrr)

# Loading all functions in package directory
devtools::load_all()

# Creating simple list of basic parameters, so produce population-level
# parameter samples that can be used to then create individual-level samples
obs_sim <- prepare_sim_data(P = 20)

# Generating population-level parameter samples in the right form for 
# then generating individual-level samples
dt_ind_params_proc <- prepare_ind_params(obs_sim)

dt_ind_params_proc$c_0 <- 50
dt_ind_params_proc$c_p_mean <- -0.3
dt_ind_params_proc$t_p_mean <- log(5)
dt_ind_params_proc$t_lod_mean <- log(25)

# Setting known effect size estimates for the three process parameters of 
# interest: c_p, t_p and t_lod
c_p_eff_man = 0.4
t_p_eff_man = 1.4
t_lod_eff_man = 0.7

# Generating individual-level parameter samples for simulated baseline group
dt_ind_params_baseline <- simulate_ind_params(
  obs = obs_sim,
  parameters = dt_ind_params_proc,
  time_range = -5:30,
  c_p_eff = 1,
  t_p_eff = 1, 
  t_lod_eff = 1
)

# Generating individual-level parameter samples for simulated adjusted group
dt_ind_params_adjusted <- simulate_ind_params(
  obs = obs_sim,
  parameters = dt_ind_params_proc,
  time_range = -5:30,
  c_p_eff = c_p_eff_man,
  t_p_eff = t_p_eff_man, 
  t_lod_eff = t_lod_eff_man, 
  covariate_description = "adjusted_group"
)

# Putting both data.tables together
dt_ind_params <- rbind(
  dt_ind_params_baseline,
  dt_ind_params_adjusted)

# Remake IDs so that they are unique after combining the data from
# the two separate groups
dt_ind_params[, id := .GRP,
              by = c("id", "covariate_1")]

# Simulating ct trajectories from the individual-level parameter samples 
dt_trajs <- simulate_cts(
  dt_ind_params, by = c("t", "id", "covariate_1"))

# Plotting the resulting datasets
dt_trajs |>
  ggplot() +
  geom_point(aes(x = t, y = ct_value, colour = covariate_1)) +
  geom_line(aes(x = t, y = ct_value, colour = covariate_1))

# Indexing time by the first positive test, rather than exposure time,
# to make sampling the time points of the trajectory and interpreting
# the samples as raw data more simple
dt_trajs <- index_by_first_positive(dt_trajs)
dt_trajs[, onset_time := as.integer(onset_time - t_first_pos)]

# Sampling a number of points from each individuals trajectory
# can specify the lower and upper bounds (of an assumed uniform distribution)
# for the number of positive and negative samples per individual
dt_data_sample <- sample_from_trajectories(
  dt_trajs, min_pos = 2, max_pos = 8, min_neg = 1, max_neg = 2)

# Redoing the levels of the covariate factor so that the baseline group is first
# as this is how the baseline is chosen within the stan code
dt_data_sample[
  , covariate_1 := forcats::fct_relevel(
    covariate_1, "baseline")]

# Setting which parameters to adjust
adj_params <- c(
  "t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")

# Specify the CT summary parameter design matrix
ct_model_sim <- subject_design(
  ~ 1 + covariate_1,
  data = dt_data_sample,
  params = adj_params,
  preds_sd = 0.2
)

# Compiling model
mod <- cmdstan_model(
  "stan/ct_trajectory_model.stan",
  include_paths = "stan",
  stanc_options = list("O1")
)

# Getting the data into the correct format for Stan
sim_stan_data <- epict_to_stan(
  dt_data_sample, 
  ct_model = ct_model_sim,
  individual_variation = 0.025,
  switch = FALSE,
  onsets = TRUE)

# Fitting the model
fit_sim <- mod$sample(
  data = sim_stan_data,
  # init = stan_inits(sim_stan_data),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Extracting fitted parameter draws
draws <- extract_draws(fit_sim)

# Adjusting them by the inferred covariate effect sizes
dt_posterior_samples_wide <- adjust_params(
  draws, design = ct_model_sim$design,
  onsets_flag = TRUE)

# Adding the baseline category to the draws
dt_posterior_samples_wide <- add_baseline_to_draws(
  dt_posterior_samples_wide,
  "baseline",
  onsets_flag = TRUE)

# setnames(dt_posterior_samples_wide,
#          c("t_p_mean", "c_p_mean", "t_lod_mean"), 
#          c("t_p", "c_p", "t_lod"))

# Transforming back to the natural scale
dt_posterior_samples_wide <- dt_posterior_samples_wide |> 
  transform_to_natural()

# Removing unused parameter names
dt_posterior_samples_wide[, t_s := NULL]
dt_posterior_samples_wide[, c_s := NULL]
dt_posterior_samples_wide[, inc_mean := NULL]
dt_posterior_samples_wide[, inc_sd := NULL]
dt_posterior_samples_wide[, .chain := NULL]
dt_posterior_samples_wide[, .iteration := NULL]

# Melting into long format for plotting  
dt_posterior_samples_long <- melt(
  dt_posterior_samples_wide,
  id.vars = c(".draw", "predictor"), 
  variable.name = "parameter")

# Renaming variables
dt_posterior_samples_long[, .draw := NULL]
dt_posterior_samples_long[predictor == "covariate_1adjusted_group",
                          predictor := "adjusted"]

# Sampling from the priors for plotting
dt_prior_samples <- sample_pop_priors(10000, scale_type = "natural")
dt_prior_samples[, t_s := NULL]
dt_prior_samples[, c_s := NULL]
dt_prior_samples[, predictor := "prior"]

# Combining prior and posterior samples
dt_samples <- rbind(dt_prior_samples, dt_posterior_samples_long)

# Building table of draws
dt_params_true_baseline <- data.table(
  t_p = dt_ind_params_proc$t_p_mean,
  c_0 = dt_ind_params_proc$c_0,
  c_p = dt_ind_params_proc$c_p_mean,
  t_lod = dt_ind_params_proc$t_lod_mean,
  sigma = dt_ind_params_proc$sigma,
  predictor = "baseline") 

# Making known parameter dataframe a data.table
dt_params_true_adj <- data.table()

# Adjusting by the assumed to be true effect sizes for each parameter
dt_params_true_adj[, c_0 := dt_params_true_baseline$c_0]
dt_params_true_adj[, c_p := dt_params_true_baseline$c_p*c_p_eff_man]
dt_params_true_adj[, t_p := dt_params_true_baseline$t_p*t_p_eff_man]
dt_params_true_adj[, t_lod := dt_params_true_baseline$t_lod*t_lod_eff_man]
dt_params_true_adj[, sigma := dt_params_true_baseline$sigma]
dt_params_true_adj[, predictor := "adjusted"]

# Combining the baseline and assumed to be true draws
dt_params_true <- rbind(
  dt_params_true_baseline, dt_params_true_adj) |> 
  transform_to_natural()

# Removing the switch parameter samples
dt_params_true[, t_s := NULL]
dt_params_true[, c_s := NULL]

# Melting into long format
dt_params_true_long <- melt(dt_params_true,
                            id.vars = "predictor", 
                            variable.name = "parameter")

# Releveling factors for plots
dt_params_true_long[, predictor := forcats::fct_relevel(predictor, "baseline")]

# Naming the parameters for the plot
params_plot <- c(
  "Ct value at peak", "Timing of peak", "Timing LOD reached")

# Renaming parameters for plot
dt_samples[parameter == "c_p", parameter := "Ct value at peak"]
dt_samples[parameter == "t_p", parameter := "Timing of peak"]
dt_samples[parameter == "t_lod", parameter := "Timing LOD reached"]

# Renaming parameters for plot
dt_params_true_long[parameter == "c_p", parameter := "Ct value at peak"]
dt_params_true_long[parameter == "t_p", parameter := "Timing of peak"]
dt_params_true_long[parameter == "t_lod", parameter := "Timing LOD reached"]

# Plotting
p_sim_no_priors <- ggplot() + 
  geom_density(data = dt_samples[predictor != "prior" & parameter %in% params_plot & value < 35],
               aes(x = value, fill = predictor)) + 
  # geom_density(data = dt_samples[predictor == "prior" & parameter %in% params_plot],
  #              aes(x = value), alpha = 0.25, colour = "grey", fill = "grey") +
  geom_vline(data = dt_params_true_long[parameter %in% params_plot],
             aes(xintercept = value, colour = predictor),
             linetype = "dashed") +
  theme_minimal() + 
  labs(title = "Without priors", tag = "A") +
  facet_wrap(~parameter, scales = "free")

p_sim_with_priors <- ggplot() + 
  geom_density(data = dt_samples[predictor != "prior" & parameter %in% params_plot & value < 35],
               aes(x = value, fill = predictor)) + 
  geom_density(data = dt_samples[predictor == "prior" & parameter %in% params_plot],
               aes(x = value), alpha = 0.25, colour = "grey", fill = "grey") +
  geom_vline(data = dt_params_true_long[parameter %in% params_plot],
             aes(xintercept = value, colour = predictor),
             linetype = "dashed") +
  theme_minimal() + 
  labs(title = "With priors", tag = "B") +
  facet_wrap(~parameter, scales = "free")

# Putting the two panels together
p_sim <- p_sim_no_priors/p_sim_with_priors +
  plot_layout(guides = "collect")

# Saving the data required to remake the figure
saveRDS(dt_samples, "outputs/plot_data/supplement/figure_S26_samples.rds")
saveRDS(dt_params_true_long, "outputs/plot_data/supplement/figure_S26_true_params.rds")

# Saving the plots
ggsave("outputs/figures/pdfs/supplement/figure_S26.pdf", 
       p_sim,
       width = 8,
       height = 5,
       bg = "white")

ggsave("outputs/figures/pngs/supplement/figure_S26.png", 
       p_sim,
       width = 8,
       height = 5,
       bg = "white")
