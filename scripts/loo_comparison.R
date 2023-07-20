library(loo)

# running the setup script
source("scripts/setup_loo_comparison.R")

# Fit the model
fit_no_voc <- epict(
  obs,
  ct_model = ct_model_no_voc,
  adjustment_model  = adjustment_model,
  likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
  onsets = TRUE, # include the symptom onset data and likelihood components
  switch = FALSE, # include the longer wane part of the model
  individual_variation = 0.025, # control the amount of 
  individual_correlation = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_loglik = TRUE,
  seed_manual = 123
)

fit_voc <- epict(
  obs,
  ct_model = ct_model_voc,
  adjustment_model  = adjustment_model,
  likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
  onsets = TRUE, # include the symptom onset data and likelihood components
  switch = FALSE, # include the longer wane part of the model
  individual_variation = 0.025, # control the amount of 
  individual_correlation = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_loglik = TRUE,
  seed_manual = 123
)

# saves fitted object - usually results in a very large file (> 100Mbs)
# if the full model is fit
# fit_no_voc$save_object("outputs/fits/fit_no_voc.rds")
# fit_voc$save_object("outputs/fits/fit_voc.rds")

# load in fits
fit_without_voc <- readRDS("outputs/fits/fit_no_voc.rds")
fit_with_voc <- readRDS("outputs/fits/fit_voc.rds")

loo_no_voc <- fit_no_voc$loo()
loo_voc <- fit_voc$loo()

loo_compare(loo_no_voc, loo_voc)

