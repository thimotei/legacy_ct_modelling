# sourcing the setup file, to load in the data, the model structure, etc
source("scripts/setup/main.R")

# seed used for final fits
seed_eff <- 485764
set.seed(seed_eff)

# Fit the model if it hasn't been run and no fit object exists locally
fit_no_onsets <- epict(
  obs,
  ct_model = ct_model,
  adjustment_model  = adjustment_model,
  likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
  onsets = FALSE, # include the symptom onset data and likelihood components
  switch = FALSE, # include the longer wane part of the model
  individual_variation = 0.025, # control the amount of
  individual_correlation = 2,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 1000,
  iter_sampling = 500,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_loglik = FALSE,
  informative_priors = TRUE
)

# saves fitted object - usually results in a very large file (> 100Mbs)
# if the full model is fit
fit_no_onsets$save_object("outputs/fits/fit_no_onsets.rds")
