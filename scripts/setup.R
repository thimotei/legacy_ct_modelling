# running inference 
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
library(ggnewscale)
library(cowplot)
library(forcats)

# loading all functions in package directory
devtools::load_all()

# load in data processed in scripts/process-data.R
dt_clean <- readRDS(here("data/processed-data.rds"))

# Do additional processing to filter for the desired number of swabs
# per positive episode
obs <- subset_data(dt_clean,
                   no_pos_swabs = 2,
                   first_pos_min = -15)

# Plot the raw data
p1_raw <- plot_obs(obs, col = factor(swab_type)) +
  labs(col = "Swab type") +
  facet_wrap(~ct_type)

# Specify which params adjusting for (see params_avail_to_adjust() for options)
# Here all available options (can also specify this using "all")
adj_params <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")

# Specify the CT summary parameter design matrix
ct_model <- subject_design(
  ~ 1 + VOC + symptoms + no_exposures + age_group + t_since_last_exposure,
  data = obs,
  params = adj_params,
  preds_sd = 0.2
)

# Specify the model to use to adjust Cts globally - we use this to adjust for
# swab type and the gene target: ORF1AB, N gene and S gene (where available)
adjustment_model <- test_design(
  ~ 1 + swab_type + ct_type,
  data = obs,
  preds_sd = 1
)

