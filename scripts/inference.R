# running the setup script
source("scripts/setup.R")

# Fit the model
fit <- epict(
  obs,
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

# saves fitted object - usually results in a very large file (> 100Mbs)
# if the full model is fit
fit$save_object("outputs/fits/fit_full.rds")
