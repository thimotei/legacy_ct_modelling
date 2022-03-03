transform_to_model <- function(draws) {
  draws <- data.table::copy(draws)

  draws[,
    `:=`(
      t_p = exp(t_p_mean),
      t_s = exp(t_s_mean),
      t_lod = exp(t_lod_mean),
      c_0 = c_0,
      c_s = c_0 * plogis(c_s_mean)
    )
  ][,
    c_p := c_s * plogis(c_p_mean)
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
  return(draws)
}