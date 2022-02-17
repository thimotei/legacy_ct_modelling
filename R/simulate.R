obs <- list(
  P = 20,
  any_onsets = 1,
  onset_time = rep(0, 20),
  c_lod = 40,
  lmean = get_inc_period()$inc_mean_p,
  lsd = get_inc_period()$inc_sd_p
)

simulate_obs <- function(obs = obs,
                         parameters = stan_inits(obs)(),
                         time_range = -10:30,
                         sample_density = 2:8) {

  params <- with(parameters,
    data.table::data.table(
      id = 1:obs$P,
      onset_time = obs$onset_time,
      T_e = T_e,
      t_p = exp(t_p_mean + t_p_var * t_p_raw),
      t_s = exp(t_s_mean + t_s_var * t_s_raw),
      t_lod = exp(t_lod_mean + t_lod_var * t_lod_raw),
      c_0 = c_0,
      c_s = c_0 * plogis(c_s_mean + c_s_var * c_s_raw)
    )[,
      c_p := c_s * plogis(c_p_mean + c_p_var * c_p_raw)
    ][,
      c_lod := obs$c_lod,
    ][,
      t_lod_abs := t_p + t_s + t_lod
    ][,
      sigma := sigma
    ]
  )

  times <- data.table::data.table(
    t = time_range
  )[,
    sample := 1:.N
  ]

  ct_trajs <- merge(
    params[, tid := 1], times[, tid := 1], by = "tid",
    allow.cartesian = TRUE
  )[,
    tid := NULL][,
    exp_ct := ct_hinge_long(
        t, c0 = c_0, cp = c_p, cs = c_s, clod = c_0, te = 0,
        tp = t_p, ts = t_s, tlod = t_lod
      ),
    by = c("id", "t", "sample")
  ][,
    ct_value := rtruncnorm(
        1, b = c_0, mean = exp_ct, sd = sigma
    )
  ][ct_value >= c_lod,
    ct_value := c_lod
  ][,
    pcr_res := ifelse(ct_value < c_lod, 1, 0)
  ]

  ct_trajs <- index_by_first_positive(ct_trajs)

  if (!is.null(sample_density)) {
    ct_trajs <- ct_trajs[,
     .SD[sample %in% sample(.N, sample(sample_density, 1))], by = "id"
    ]
  }

  return(ct_trajs)
}
