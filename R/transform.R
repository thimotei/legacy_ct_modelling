transform_to_model <- function(draws) {
  draws <- data.table::copy(draws)

  draws[,
    `:=`(
      t_p = exp(t_p),
      t_s = exp(t_s),
      t_lod = exp(t_lod),
      c_0 = c_0,
      c_s = c_0 * plogis(c_s)
    )
  ][,
    c_p := c_s * plogis(c_p)
  ]
  return(draws[])
}

transform_to_natural <- function(draws) {
  draws <- transform_to_model(draws)
  draws[,
    t_lod :=  t_p + t_s + t_lod
  ][,
    t_s := t_p + t_s,
  ]
  return(draws[])
}

transform_ip_to_natural <- function(draws) {
  draws <- draws[, `:=`(
    nat_inc_mean = exp(inc_mean + (inc_sd^2) / 2),
    nat_inc_sd = sqrt((exp(inc_sd^2) - 1) * exp(2 * inc_mean + inc_sd^2))
  )]
  return(draws[])
}