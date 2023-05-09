#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(stringr)
library(purrr)
library(cowplot)

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
  preds = 0,
  ct_preds = 0,
  K = 5,
  switch = 1,
  ind_var_sd = 0.025
)

# Simulate from the centre of the prior for all parameters
# based on initial conditions used for the stan model.
params_sample_1 <- simulate_params(
  obs = obs,
  parameters = epict_inits(obs)(),
  time_range = -1:25,
  t_p_eff_size = 1.0)

ct_sample_1 <- simulate_obs(params_sample_1, sample_density = 4:8)

ct_sample_1[,VOC := "Delta"][, 
  swab_type := "Dry"]

# Simulate from the centre of the prior for all parameters
# based on initial conditions used for the stan model.
params_sample_2 <- simulate_params(
  obs = obs,
  parameters = epict_inits(obs)(),
  time_range = -1:25,
  t_p_eff_size = 1.001)

ct_sample_2 <- simulate_obs(params_sample_2, sample_density = 4:8)

ct_sample_2[,VOC := "Omicron (BA.1)"][, 
  swab_type := "Wet"][, id := id + obs$P]

ct_sample <- rbindlist(list(ct_sample_1, ct_sample_2))[,
  VOC := factor(VOC)][, swab_type := factor(swab_type)]

# Specify which params adjusting for (see params_avail_to_adjust() for options)
# Here all available options (can also specify this using "all")
adj_params <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")

# Specify the CT summary parameter design matrix
ct_model <- subject_design(
  ~ 1 + VOC,
  data = ct_sample,
  params = adj_params,
  preds_sd = 0.2
)

# Specify the model to use to adjust Cts globally - we use this to adjust for
# swab type and the gene target: ORF1AB, N gene and S gene (where available)
adjustment_model <- test_design(
  ~ 1 + swab_type,
  data = ct_sample,
  preds_sd = 1
)

# plot of subset of data
plot_obs(ct_sample) +
  facet_wrap(vars(factor(id), VOC))

# Fit the model
fit <- epict(
  ct_sample,
  ct_model = ct_model,
  adjustment_model  = adjustment_model,
  likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
  onsets = TRUE, # include the symptom onset data and likelihood components
  switch = FALSE, # include the longer wane part of the model
  individual_variation = 0.025, # control the amount of 
  individual_correlation = 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_loglik = FALSE
)

draws <- extract_draws(fit)

adj_draws <- adjust_params(draws, 
                           design = ct_model$design, 
                           onsets_flag = TRUE) 


adj_draws <- add_baseline_to_draws(adj_draws, 
                                   "Delta", 
                                   onsets_flag = TRUE) %>% 
  transform_to_natural()

adj_draws_long <- melt(adj_draws[, c(".draw", "c_0", "c_p", 
                                     "t_p", "t_lod",
                                     "inc_mean", "inc_sd",
                                     "predictor")],
                       id.vars = c(".draw", "predictor"),
                       variable.name = "parameter") 

ind_level_params <- c("c_0", "c_p", "c_s", "c_lod", "t_p")

adj_draws_long %>% 
  ggplot() +
  geom_density(aes(x = value, fill = predictor), alpha = 0.2) +
  geom_vline(data = ct_sample_nat, 
             aes(xintercept = mean(c_p)), linetype = "dashed") +
  facet_wrap(~parameter, scales = "free")

ct_sample_nat_pop <- ct_sample_nat[, c("t_p_mean", "t_s_mean", 
                                       "t_lod_mean", "c_p_mean", "c_s_mean")]

pop_params <- ("c_p_mean", "t_p_mean")

ct_sample_nat <- ct_sample %>% 
  transform_to_natural()

# # compiling model
# mod <- cmdstan_model(
#   "stan/ct_trajectory_model.stan",
#   include_paths = "stan",
#   stanc_options = list("O1")
# )
# 
# sim_stan_data <- epict_to_stan(ct_sample, onset = TRUE)
# 
# # fitting the model - not very quick, as many iterations hit the
# # max_tree_depth at the moment
# fit_sim <- mod$sample(
#   data = sim_stan_data,
#   init = stan_inits(sim_stan_data),
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 2000
# )

# Extract and plot posterior predictions
sim_pp_plot <- pp_plot <- plot_obs(
  obs = ct_sample,
  ct_traj =  extract_ct_trajectories(fit_sim),
  pp = summarise_pp(fit_sim, ct_sample),
  samples = 10, traj_alpha = 0.05,
  col = factor(swab_type)
) +
  labs(col = "Swab type")
  facet_wrap(vars(factor(id)))
  facet_wrap(vars(factor(id))) +

ggsave("outputs/figures/sim_pp.png", sim_pp_plot, height = 10, width = 10)

# Extract and plot population level posterior predictions for the CT model
sim_draws <- extract_draws(fit_sim)

adj_draws <- adjust_draws(sim_draws)

sim_parameter_pp <- plot_summary(sim_draws)

ggsave(
  "outputs/figures/sim_parameter_pp.png",
  sim_parameter_pp, width = 12, height = 8,
)


ggsave(
  "outputs/figures/sim_population_ct_pp.png",
  sim_pop_pp, width = 8, height = 8,
)
