summarise_pop_pp <- function(fit) {
  draws <- fit$summary(
    variables = c(
      "t_p_mean", "t_s_mean", "t_lod_mean", "c_p_mean",
      "c_s_mean", "inc_mean", "inc_sd"
    )
  )
  draws <- data.table::as.data.table(draws)
  return(draws[])
}

summarise_coeff_pp <- function(fit, params, exponentiate = FALSE) {

  params <- params_avail_to_adjust(params)
  params <- names(params[purrr::map_lgl(params, ~ . == 1)])

  draws <- fit$summary(
    variables = paste0("beta_", params)
  )
  draws <- data.table::as.data.table(draws)

  if (exponentiate) {
    beta_cols <-  c("mean", "median", "sd", "mad", "q5", "q95")
    draws[, (beta_cols) := lapply(.SD, exp), .SDcols = beta_cols]
  }
  return(draws[])
}

summarise_effects <- function(draws, design, variables, exponentiate = TRUE) {
  eff_draws <- extract_coeffs(
      draws, exponentiate = exponentiate, design = design, variables
    )

    by <- "variable"
    if (!missing(design)) {
      by <- c(by, "predictor")
    }
    eff_summary <- summarise_draws(eff_draws, by = by)
    return(eff_summary)
}

summarise_draws <- function(draws, by = c("id", "time")) {

  out <- draws[,
      .(median = quantile(value, c(0.5), na.rm = TRUE),
        lo90 = quantile(value, c(0.05), na.rm = TRUE),
        lo60 = quantile(value, c(0.20), na.rm = TRUE),
        lo30 = quantile(value, c(0.35), na.rm = TRUE),
        hi30 = quantile(value, c(0.65), na.rm = TRUE),
        hi60 = quantile(value, c(0.80), na.rm = TRUE),
        hi90 = quantile(value, c(0.95), na.rm = TRUE)
      ),
      by = by
    ]
  return(out[])
}
