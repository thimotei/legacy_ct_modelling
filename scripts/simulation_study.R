#simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(cowplot)
library(stringr)

# loading all functions in package directory
devtools::load_all()

# example of a single trajectory
t_max <- 30
t_step <- 1
t_input <- seq(1, t_max, t_step)

# free parameters in the model cp, cs, te, tp, ts, tlod
# fixed parameters in mechanistic model are c0 and clod, both fixed at 40
c0 <- 40
clod <- 40

# fixed parameters in statistical model are sigma_obs
sigma_obs_known <- 1

n <- 10 # total number of individuals being simulated

# simulating trajectories. all parameters are sampled from uniform distributions
# where the minimum and maximum of each are arguments of the simulating
# function
ext_ct_dt <- simulate_ct_trajectories(t_max = 30, t_stepsize = 1, 
                                      cp_min = 0.3, cp_max = 0.7,
                                      cs_min = 0.2, cs_max = 0.8,
                                      te_min = 1, te_max = 7, 
                                      tp_min = 1, tp_max =7,
                                      ts_min = 1, ts_max = 7, 
                                      tlod_min = 5, tlod_max = 10,
                                      sigma_obs = 1)

# quick plot of simulated data
ext_ct_dt %>% 
  ggplot() + 
  geom_point(aes(x = t, y = ct_value_std, colour = pcr_res)) +
  facet_wrap(vars(id)) + 
  custom_plot_theme()

# saving the true parameters for comparison to estimated values
true_params <- ext_ct_dt[, .(id = unique(id), 
                             cs = unique(cs),
                             cp = unique(cp),
                             te = unique(te),
                             tp = unique(tp),
                             ts = unique(ts),
                             tlod = unique(tlod))]

# sampling a "realistic size" subset of the data. I.e. between 3 and 8 samples 
# at random times per person
ext_ct_dt_sample <- ext_ct_dt[, .SD[t %in% sample(.N, sample(3:8, 1))],
                              by = "id"]

# quick plot of subset of data
ext_ct_dt_sample %>% 
  ggplot(aes(x = t, y = exp_ct_std)) + 
  geom_point() +
  facet_wrap(vars(id)) + 
  custom_plot_theme()

# compiling model to test inference
mod <- cmdstan_model("stan/ct_trajectory_model_individual.stan",
                     include_paths = "~/lshtm/legacy_ct_modelling/stan")

#--- running inference
n.chains <- 4
stan_data_simulated <- stan_data_fun(ext_ct_dt, simulated = TRUE)
options(mc.cores = 8)

# fitting the model - not very quick, as many iterations hit the
# max_tree_depth at the moment
fit_sim <- mod$sample(
  data = stan_data_simulated,
  seed = 123,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  init = 0.1
)

# extracting draws and putting them nicely into a data.table
draws_dt <- as.data.table(fit_sim$draws())

# extracting Ct fits. Bit slow as it is at the moment
ct_dt_draws <- extract_ct_fits(draws_dt[variable %like% "ct"])

# summarising trajectories using median and 95% CrI
ct_dt_draws_summary <- ct_trajectory_summarise(ct_dt_draws)

# plotting summaries of fitted trajectories against simulated data
plot_ct_trajectories(ct_dt_draws_summary, ext_ct_dt)


