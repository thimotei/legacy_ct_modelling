update_variable_labels <- function(draws, reverse = FALSE) {

draws <- data.table::copy(draws)
params <- c(
  "c_0", "c_p", "c_s","t_p", "t_s", "t_lod", "inc_mean", "inc_sd",
  "nat_inc_mean", "nat_inc_sd", "ct_shift", "ct_scale"
)

clean_params <- c(
  "Ct value at latent limit of detection",
  "Ct value at peak",
  "Ct value at switch",
  "Time of peak",
  "Time of switch",
  "Time of latent limit of detection",
  "Incubation period (log mean)",
  "Incubation period (log sd)",
  "Incubation period (mean)",
  "Incubation period (sd)",
  "Ct intercept adjustment",
  "Ct multiplicative adjustment"
)

if (reverse) {
  params <- rev(params)
  clean_params <- rev(clean_params)
}

  draws <- draws[variable %in% params][,
    variable := factor(
      variable,
      levels = params,
      labels = clean_params
    )
  ]
  return(draws[])
}