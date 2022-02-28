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
library(tidyr)

# loading all functions in package directory
files <- list.files("R", "*.R", full.names = TRUE)
walk(files, source)

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

#--- Drop patients with gaps between positive tests of more than 60 days
ids_spurious_gaps <- dt_clean[,
  .(spurious = any(abs(time_since_first_pos) > 60) || any(abs(onset_time) > 60)
   ), by = "id"
]
ids_spurious_gaps[spurious == TRUE][]
dt_clean <- dt_clean[ids_spurious_gaps, on = "id"]
dt_clean <- dt_clean[spurious == FALSE | is.na(spurious)]

#--- Drop subject with ID 978 due to large mismatch between onset and positive
#--  tests
dt_clean <- dt_clean[id != 978]

#--- quick plots of the raw data, stratified by variant
dt_2_tests <- subset_data(
  dt_clean, voc = c("Delta", "Omicron"), no_pos_swabs = 2
)
dt_delta <- subset_data(dt_clean, voc = c("Delta"), no_pos_swabs = 2)
dt_omicron <- subset_data(dt_clean, voc = "Omicron", no_pos_swabs = 2)

p1_raw <- plot_obs_ct(dt_2_tests)

#--- compile models
mod <- cmdstan_model("stan/ct_trajectory_model.stan", include_paths = "stan")

#--- fit
stan_data <- data_to_stan(dt_2_tests, likelihood = TRUE, onsets = TRUE)
stan_data_delta <- data_to_stan(dt_delta, likelihood = TRUE, onsets = TRUE)
stan_data_omicron <- data_to_stan(dt_omicron, likelihood = TRUE, onsets = TRUE)

fit <- mod$sample(
  data = stan_data,
  init = stan_inits(stan_data),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

fit_delta <- mod$sample(
  data = stan_data_delta,
  init = stan_inits(stan_data_delta),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)


fit_omicron <- mod$sample(
  data = stan_data_omicron,
  init = stan_inits(stan_data_omicron),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

#--- extract and plot posterior predictions
ct_draws <- extract_ct_trajectories(fit)

ct_summary <- summarise_draws(
  copy(ct_draws)[,
    time_since_first_pos := as.integer(time_since_first_pos)
    ],
  by = c("id", "time_since_first_pos")
)

ct_pp <- extract_posterior_predictions(fit, dt_2_tests)
ct_pp <- summarise_draws(
  ct_pp[, value := sim_ct], by = c("id", "t", "pcr_res", "obs")
)

# plotting summaries of fitted trajectories against simulated data
pp_plot <- plot_obs_ct(
  dt_2_tests, ct_draws[iteration <= 10], ct_pp, traj_alpha = 0.05
)

ggsave("outputs/figures/pp.png", pp_plot, height = 16, width = 16)

#--- extract and plot population level posterior predictions

pop_draws <- extract_draws(fit)

pop_ct_draws <- transform_to_model(pop_draws) %>%
  simulate_cts(time_range = 0:60, obs_noise = FALSE)


# Extract population level CT parameter samples

# Make population level expected CT trajectories

# Plot expected CT trajectories

# Plot population level CT parameter samples

# Joint plot of expected CT trajectories and population level CT parameter samples



draws <- rbind(as.data.table(fit_delta$draws())[, voc := "delta"],
               as.data.table(fit_omicron$draws())[, voc := "omicron"])


#--- gathering and plotting population-level posteriors by VOC
pop_params_rel <- c("c_0", "c_p_mean", "c_s_mean", 
                    "t_p_mean", "t_s_mean", "t_lod_mean")

pop_posteriors_wide <- transform_pop_posteriors(pop_params_rel, draws)

pop_params_abs <- c("c_0_abs", "c_p_mean_abs", "c_s_mean_abs", 
                    "t_p_mean_abs", "t_s_mean_abs", "t_lod_mean_abs")

pop_params_abs_labs <- c("Ct value at zero & LOD", "Ct value at peak", 
                         "Ct value at switch", "Time of peak", "Time of switch", 
                         "Time of limit of detection")

pop_params_abs_melt <- c("iteration", "voc", pop_params_abs)

pop_posteriors_long <- melt(pop_posteriors_wide[, ..pop_params_abs_melt], 
                            measure.vars = pop_params_abs)

p_pop_posteriors <- plot_pop_posteriors(
  pop_posteriors_long %>% 
    setnames(., "voc", "VOC") %>% 
    .[, variable := factor(variable, 
                           levels = pop_params_abs,
                           labels = pop_params_abs_labs)]
  )


#--- gathering and plotting population-level posterior predictive curves by VOC
pop_pp_samples_dt <- pop_pp_samples(pop_posteriors_wide, t_max = 40, t_step = 0.1)
pop_pp_summary_dt <- summarise_draws(pop_pp_samples_dt, by = c("time", "voc"))

p_pred_summary <- plot_pop_pp(pop_pp_summary_dt,
                              pop_pp_samples_dt,
                              no_samples = 100) 


#--- combining posteriors and posterior predictive plots
p_final <- p_pop_posteriors / p_pred_summary + plot_layout(guides = "collect")
  
ggsave("outputs/figures/posteriors_and_predictive.png",
       p_final,
       width = 8,
       height = 8,
       bg = "white")
