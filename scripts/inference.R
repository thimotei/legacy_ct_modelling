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
obs <- subset_data(dt_clean, no_pos_swabs = 2)

# Plot the raw data
p1_raw <- plot_obs(obs, col = factor(swab_type)) +
  labs(col = "Swab type")

# Specify which params adjusting for (see params_avail_to_adjust() for options)
# Here all available options (can also specify this using "all")
adj_params <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")

# Specify the CT summary parameter design matrix
ct_model <- subject_design(
  ~ 1 + VOC + symptoms + no_vaccines + time_since_last_dose,
  data = obs,
  params = adj_params,
  preds_sd = 0.4
)

# Specify the model to use to adjust CTs globally
adjustment_model <- test_design(
  ~ 1 + swab_type,
  data = obs,
  preds_sd = 1
)

# Add a function mapping labels to each term
update_predictor_labels <- function(dt) {
  dt <- dt[,
    predictor := factor(
      predictor,
      levels =  c(
        "no_vaccines2", "symptomsasymptomatic", "VOCDelta", "VOCBA.2",
        "swab_typeVTM", "time_since_last_dose"
      ),
      labels = c(
        "2 vaccines", "asymptomatic", "Delta", "BA.2",
        "Swab type", "Time since last dose"
      )
    )
  ]
  return(dt[])
}

# Fit the model
fit <- epict(
  obs,
  ct_model = ct_model,
  adjustment_model  = adjustment_model,
  likelihood = TRUE,
  onsets = TRUE,
  switch = FALSE,
  individual_variation = 0.2,
  individual_correlation = 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12
)

# Extract and plot posterior predictions
pp_plot <- plot_obs(
  obs = obs,
  pp = summarise_pp(fit, obs),
  samples = 10, traj_alpha = 0.05,
  col = factor(swab_type)
) +
  labs(col = "Swab type") +
  facet_wrap(vars(factor(id)))

ggsave(
  "outputs/figures/pp.png", pp_plot, height = 16, width = 16
)


# Extract posterior predictions
draws <- extract_draws(fit)

# Extract effect sizes and make a summary plot
eff_plot <- plot_effect_summary(
  draws, ct_design = ct_model$design, variables = adj_params,
  adjustment_design = adjustment_model$design,
  variable_labels = update_variable_labels,
  col = predictor, position = position_dodge(width = 0.6)
) &
  scale_colour_brewer(palette = "Dark2") &
  labs(col = "Adjustment")

ggsave(
  "outputs/figures/effects.png",
  eff_plot, width = 9, height = 12,
)

# Add adjusted effects to draws
adj_draws <- adjust_params(draws, design = ct_model$design)

# Filter for just adjustments that summary shows appear to differ from base case
adj_draws <- adj_draws[
  predictor %in% c("no_vaccines2", "symptomsasymptomatic",  "VOCDelta")
] %>%
  update_predictor_labels()

adj_draws <- rbind(
  extract_param_draws(draws)[,
   predictor := "Baseline (Omicron & symptomatic & 3 vaccines)"
  ],
  adj_draws
)[, predictor := factor(predictor)]

# Extract and plot population level posterior predictions for the CT model
# for the baseline case and each adjusted case individually
parameter_pp <- plot_summary(
  adj_draws, fill = predictor, colour = predictor, by = "predictor",
  simulated_samples = 1000, samples = 0,
  ct_time_range = seq(0, 30, by = 0.25), ip_time_range = seq(0, 15, by = 0.25)
) &
  scale_colour_brewer(palette = "Dark2") &
  scale_fill_brewer(palette = "Dark2") &
  labs(fill = "Adjustment", colour = "Adjustment")

ggsave(
  "outputs/figures/parameter_pp.png",
  parameter_pp, width = 12, height = 8,
)
