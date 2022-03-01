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

summarise_draws <- function(draws, by = c("id", "time")) {

  out <- draws[,
      .(median = quantile(value, c(0.5)),
        lo90 = quantile(value, c(0.05)),
        lo60 = quantile(value, c(0.20)),
        hi60 = quantile(value, c(0.80)),
        hi90 = quantile(value, c(0.95))
      ),
      by = by
    ]
  return(out[])
}
