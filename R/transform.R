transform_to_model <- function(draws, onsets_flag = FALSE) {
  draws <- data.table::copy(draws)

  if (is.null(draws[["t_s"]])) {
    draws[, t_s := -Inf]
  }

  if (is.null(draws[["c_s"]])) {
    draws[, c_s := c_0]
  }
  if(onsets_flag == TRUE) {
    draws[,
      `:=`(
        t_p = exp(t_p),
        t_s = exp(t_s),
        t_lod = exp(t_lod),
        c_0 = c_0,
        c_s = c_0 * plogis(c_s),
        inc_mean_nat = exp(inc_mean),
        inc_sd_nat = exp(inc_sd)
      )
    ][,
      c_p := c_s * plogis(c_p)
    ]
  }
  else if(onsets_flag == FALSE) {
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
    
  }
  return(draws[])
}

transform_to_natural <- function(draws) {
  draws <- transform_to_model(draws)
  draws[,
    t_lod := t_p + t_s + t_lod
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

ct_to_vl <- function(ct_value, alpha, kappa, ml_conversion){
  
  vl_out <- exp((40 - ct_value - kappa)/alpha)*ml_conversion
  return(vl_out)
}

