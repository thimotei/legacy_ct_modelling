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
adjustment.fun <- function(alpha, beta, x) {alpha + beta*x}
adjustment.dt <- fread("data/adjustment_params.csv")

dt_raw <- fread("data/ct_values_clean.csv")

# processing data, adding machine readable dates, moving dates
# back for certain swabs, based on data collection advice and adjusting
# all wet swabs upwards according to the fitted dry vs wet adjustment
# factor
dt_clean <- process_data(dt_raw)

#--- TEMPORARY, for now, as "invalid" and "inconclusive" both have many
#--- observations consistent with positive Ct values, keeping them for now
#--- will discuss with Crick partners about this
dt_clean[result == "Invalid", result := "Positive"]
dt_clean[result == "Inconclusive", result := "Positive"]
dt_clean[result == "Negative", pcr_res := 0]
dt_clean[result == "Positive", pcr_res := 1]

# calculating minimum and maximum Ct values and rescaling
# all Ct values on a standardised scale
mn <- 0
mx <- 40
dt_clean[, ct_value_std := ct_value/mx]

#--- quick plots of the raw data, stratified by variant
dt_delta <- subset_data_fun(dt_clean, voc = "Delta", no_pos_swabs = 2)
dt_omicron <- subset_data_fun(dt_clean, voc = "Omicron", no_pos_swabs = 2)

p1_raw <- plot_empiricial_data_ind(dt_delta, dt_omicron)

options(mc.cores = parallel::detectCores()) 
#--- choose here whether you want the individual-level fits 
#--- or a pooled fit
mod <- cmdstanr::stan_model("stan/ct_trajectory_model_individual.stan")

#--- fitting for delta
stan_data_delta <- stan_data_fun(dt_delta)

fit_delta <- sampling(mod,
                      chains = 4,
                      data = stan_data_delta,
                      iter = 2000)

# extracting draws and putting them nicely into a data.table
draws_delta_dt <- as.data.table(fit_delta$draws())

# extracting Ct fits. Bit slow as it is at the moment
ct_dt_delta_draws <- extract_ct_fits(draws_delta_dt)

# round time first positive to the nearest day
ct_dt_delta_draws <- ct_dt_delta_draws[,
                           time_since_first_pos := as.integer(time_since_first_pos)
]

# summarising trajectories using median and 95% CrI
ct_dt_draws_delta_summary <- ct_trajectory_summarise(
  ct_dt_delta_draws, by = c("id", "time_since_first_pos")
)

# plotting summaries of fitted trajectories against simulated data
plot_ct_trajectories(ct_dt_delta_draws_summary, dt_delta)

