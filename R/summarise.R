summarise_pop_pp <- function(fit) {
  fit$summary(
    variables = c(
      "t_p_mean", "t_s_mean", "t_lod_mean", "c_p_mean",
      "c_s_mean", "inc_mean", "inc_sd"
    )
  )
}

summarise_coeff_pp <- function(fit, params) {

  params <- params_avail_to_adjust(params)
  params <- names(params[purrr::map_lgl(params, ~ . == 1)])

  fit$summary(
    variables = paste0("beta_", params)
  )
}