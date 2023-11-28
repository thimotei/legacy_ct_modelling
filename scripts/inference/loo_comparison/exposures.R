# running the setup script
source("scripts/setup/exposures.R")

# sampling statement to produce manual seed
# seed_manual <- sample(1:2^20, 1)

# seed used for final fits
seed_eff <- 485764
set.seed(seed_eff)

# Fit the model
fit_exposures <- epict(
  obs,
  ct_model = ct_model_exposures,
  adjustment_model = adjustment_model,
  likelihood = TRUE, # flag switching likelihood component of model on or off. Off means priors alone are sampled from
  onsets = TRUE, # include the symptom onset data and likelihood components
  switch = FALSE, # include the longer wane part of the model
  individual_variation = 0.025, # control the amount of 
  individual_correlation = 1.0,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_loglik = TRUE,
  seed_manual = seed_eff,
  informative_priors = TRUE
)

# saves fitted object - usually results in a large file (> 100Mbs)
# if the full model is fit
fit_exposures$save_object("outputs/fits/exposures.rds")
