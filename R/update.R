update_ct_variables <- function(draws, reverse = FALSE) {

draws <- data.table::copy(draws)
params <- c(
  "c_0", "c_p", "c_s","t_p", "t_s", "t_lod", "inc_mean", "inc_sd"
)

clean_params <- c(
  "Ct value at limit of detection",
  "Ct value at peak",
  "Ct value at switch",
  "Time of peak",
  "Time of switch",
  "Time of limit of detection",
  "Incubation period (mean)",
  "Incubation period (sd)"
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