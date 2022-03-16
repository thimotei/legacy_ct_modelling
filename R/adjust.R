adjust_params <- function(draws, design) {
  param_draws <- extract_param_draws(draws)
  param_draws <- melt_draws(param_draws)
  eff_draws <- extract_coeffs(
    draws, design = design, variables = params_avail_to_adjust()
  )
  data.table::setnames(eff_draws, "value", "mod")

  linked_draws <- merge(
    param_draws,
    eff_draws,
    by = c(".chain", ".iteration", ".draw", "variable"),
    all.x = TRUE
  )

  # combine effects
  linked_draws[is.na(mod), mod := 0]
  linked_draws[, value := value + mod]

  linked_draws <- dcast(
    linked_draws, .chain + .iteration + .draw + predictor ~ variable
  )

  # account for unadjusted parameters for a given predictor
  cols <- as.character(unique(param_draws$variable))

  linked_draws[, (cols) := purrr::map(
    .SD, ~ ifelse(is.na(.), .[is.na(predictor)], .)
    ),
    by = c(".chain", ".iteration", ".draw"), .SDcols = cols
  ]
  linked_draws <- linked_draws[!is.na(predictor)]
  return(linked_draws[])
}
