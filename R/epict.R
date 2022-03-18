epict <- function(obs,
                  model = load_epict_model(),
                  as_data_list = epict_to_stan,
                  inits = epict_inits,
                  ct_model = subject_design(~ 1, obs),
                  adjustment_model = test_design(~ 1, obs),
                  individual_variation = 0.2, individual_correlation = 1,
                  censoring_threshold = 40, switch = TRUE,
                  onsets = TRUE, incubation_period = get_inc_period(),
                  likelihood = TRUE, output_loglik = FALSE, ...) {
  stan_data <- as_data_list(
    obs,
    ct_model = ct_model,
    adjustment_model = adjustment_model,
    individual_variation = individual_variation,
    individual_correlation = individual_correlation,
    censoring_threshold = censoring_threshold,
    switch = switch,
    onsets = onsets, incubation_period = incubation_period,
    likelihood = likelihood, output_loglik = output_loglik
  )

  fit <- model$sample(
    data = stan_data,
    init = inits(stan_data),
    ...
  )
  return(fit)
}