#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(ggplot2)
library(truncnorm)
library(cmdstanr)
library(cowplot)
library(stringr)
library(purrr)
library(lubridate)
library(patchwork)

# loading all functions in package directory
files <- list.files("R", "*.R", full.names = TRUE)
walk(files, source)

#--- loading in adjustment factors and defining function for 
#--- dry vs wet swab adjustment. We adjust the VTM swabs downwards
adj_params <- fread("data/adjustment_params.csv")

dt_raw <- fread("data/ct_values_clean.csv")

# processing data, adding machine readable dates, moving dates
# back for certain swabs, based on data collection advice and adjusting
# all wet swabs upwards according to the fitted dry vs wet adjustment
# factor
dt_clean <- process_data(dt_raw, adj_params)

#--- TEMPORARY, for now, as "invalid" and "inconclusive" both have many
#--- observations consistent with positive Ct values, keeping them for now
#--- will discuss with Crick partners about this
dt_clean[result == "Invalid", result := "Positive"]
dt_clean[result == "Inconclusive", result := "Positive"]
dt_clean[result == "Negative", pcr_res := 0]
dt_clean[result == "Positive", pcr_res := 1]


#--- quick plots of the raw data, stratified by variant
dt_2_tests <- subset_data(
  dt_clean, voc = c("Delta", "Omicron"), no_pos_swabs = 2
)
dt_delta <- subset_data(dt_clean, voc = c("Delta"), no_pos_swabs = 2)
dt_omicron <- subset_data(dt_clean, voc = "Omicron", no_pos_swabs = 2)

p1_raw <- plot_empiricial_data_ind(dt_delta, dt_omicron)

#--- compile models
mod <- cmdstan_model("stan/ct_trajectory_model.stan", include_paths = "stan")

#--- fit
stan_data <- data_to_stan(dt_2_tests)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# extracting Ct fits. Bit slow as it is at the moment
ct_draws <- extract_ct_trajectories(fit)

# summarising trajectories using median and 95% CrI
ct_summary <- summarise_draws(
  copy(ct_draws)[,
    time_since_first_pos := as.integer(time_since_first_pos)
    ],
  by = c("id", "time_since_first_pos")
)

# plotting summaries of fitted trajectories against simulated data
pp_plot <- plot_obs_ct(dt_2_tests, ct_draws[iteration <= 10], traj_alpha = 0.05)
ggsave("outputs/figures/pp.png", pp_plot, height = 10, width = 10)