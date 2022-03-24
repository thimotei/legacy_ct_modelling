# performing model selection on width of effect size prior
effect_prior_sd <- seq(0.05, 0.5, 0.05)

# varying model design matrix by the vector of assumed sd's
ct_model <- effect_prior_sd %>% 
  purrr::map(
    ~subject_design(
      ~ 1 + VOC + symptoms + no_vaccines + time_since_last_dose,
      data = obs,
      params = adj_params,
      preds_sd = .x
    )
  )

# fitting the model 10 times with each of the parameter values
all_fits <- ct_model %>% 
  purrr::map(~epict(
    obs,
    ct_model = .,
    adjustment_model = adjustment_model,
    likelihood = TRUE,
    onsets = TRUE,
    switch = FALSE,
    individual_variation = 0.2,
    individual_correlation = 1,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    output_loglik = TRUE
  )
)

# function to return loo for a cmdstan stan_fit object
# not sure how to use the standard cmdstanr syntax to do this within
# an lapply function
compute_loo <- function(stan_fit) {
  return(stan_fit$loo(save_psis = TRUE))
}

# looping over all the fit objects and computing the loo estimates
all_loo_outputs <- lapply(all_fits, compute_loo)

# comparing all 10 loo objects at once and returning the comparison, sorted
# by best fitting model first
loo::loo_compare(all_loo_outputs)