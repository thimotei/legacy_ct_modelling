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

# Simulate trajectories for 10 individuals with a test per day
sim_ct <- simulate_ct_trajectories(
  t_max = 30, t_stepsize = 1, cp_min = 10, cp_max = 20,
  cs_min = 20, cs_max = 30, te_min = 1, te_max = 7,
  tp_min = 1, tp_max = 7, ts_min = 1, ts_max = 7,
  tlod_min = 5, tlod_max = 10, c0 = 40, clod = 40, n = 10,
  sigma_obs = 1
)

plot_obs_ct(sim_ct)

# saving the true parameters for comparison to estimated values
true_params <- sim_ct[,
  .(id = unique(id), cs = unique(cs), cp = unique(cp),
    te = unique(te), tp = unique(tp), ts = unique(ts),
    tlod = unique(tlod)
  )
]

# sampling a "realistic size" subset of the data. I.e. between 3 and 8 samples
# at random times per person
ct_sample <- sim_ct[, .SD[t %in% sample(.N, sample(3:8, 1))], by = "id"]

# get time for first positive test per person
ct_sample <- index_by_first_positive(ct_sample)

# quick plot of subset of data
plot_obs_ct(ct_sample)

# compiling model
mod <- cmdstan_model("stan/ct_trajectory_model.stan", include_paths = "stan")

sim_stan_data <- data_to_stan(ct_sample, likelihood = FALSE)

# fitting the model - not very quick, as many iterations hit the
# max_tree_depth at the moment
fit_sim <- mod$sample(
  data = sim_stan_data,
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

# plotting summaries of fitted trajectories against simulated data
sim_pp_plot <- plot_obs_ct(
  ct_sample, ct_draws[iteration <= 10], traj_alpha = 0.05
)
ggsave("outputs/figures/sim_pp.png", sim_pp_plot, height = 10, width = 10)