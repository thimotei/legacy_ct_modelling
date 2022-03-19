# Ct trajectory functions used to simulate data - the same functions
# as those in the Stan model

# for Ct trajectories with the long wane
piecewise_ct <- function(t, c0, cp, cs, clod, te, tp, ts, tlod) {
  if (ts == 0) {
    cs <- cp
  }
  if (t <= te) {
    y <- c0
  } else if (t > te && t <= te + tp) {
    y <- ((t - te) * (cp - c0)) / tp  + c0
  } else if (t > te + tp && t <= te + tp + ts) {
    y <- ((t - te - tp) * (cs - cp)) / ts  + cp
  } else if (t > te + tp + ts && t <= te + tp + ts + tlod) {
    y <- ((t - te - tp - ts) * (clod - cs)) / tlod + cs
  } else if (t > tlod) {
    y <- clod
  }
  return(y)
}

simulate_cts <- function(params, time_range = 0:30, obs_noise = TRUE) {

  if (is.null(params[["id"]])) {
    params[, id := 1:.N]
  }

  if (is.null(params[["t_s"]])) {
    params[, t_s := 0]
  }

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
    exp_ct := piecewise_ct(
        t, c0 = c_0, cp = c_p, cs = c_s, clod = c_0, te = 0,
        tp = t_p, ts = t_s, tlod = t_lod
      ),
    by = c("id", "t", "sample")
  ]

  if (obs_noise) {
    ct_trajs[,
    ct_value := truncnorm::rtruncnorm(
        1, b = c_0, mean = exp_ct, sd = sigma
    )
    ][ct_value >= c_lod,
      ct_value := c_lod
    ][,
      uncensored := ifelse(ct_value < c_lod, 1, 0)
    ]
  }else{
    ct_trajs[, ct_value := exp_ct]
  }

  return(ct_trajs[])
}

simulate_ips <- function(params, time_range = 0:30) {

  if (is.null(params[["id"]])) {
    params[, id := 1:.N]
  }

  times <- data.table::data.table(
    t = time_range
  )[,
    sample := 1:.N
  ]

  ip <- merge(
    params[, tid := 1], times[, tid := 1], by = "tid",
    allow.cartesian = TRUE
  )

  ip <- ip[,
    value := dlnorm(time_range, inc_mean, inc_sd),
    by = c("id")
  ]

  return(ip[])
}

simulate_obs <- function(obs = obs,
                         parameters = stan_inits(obs)(),
                         time_range = 0:30,
                         sample_density = 2:8) {

  params <- with(parameters,
    data.table::data.table(
      id = 1:obs$P,
      swab_type = "Dry",
      onset_time = rlnorm(obs$P, inc_mean, inc_sd),
      t_inf = t_inf,
      t_p = exp(t_p_mean + ind_var[1] * ind_eta[1, ]),
      t_s = exp(t_s_mean + ind_var[2] * ind_eta[4, ]),
      t_lod = exp(t_lod_mean + ind_var[3] * ind_eta[2, ]),
      c_0 = c_0,
      c_s = c_0 * plogis(c_s_mean + ind_var[4] * ind_eta[5, ])
    )[,
      c_p := c_s * plogis(c_p_mean + ind_var[5] * ind_eta[3, ])
    ][,
      c_lod := obs$c_lod,
    ][,
      t_lod_abs := t_p + t_s + t_lod
    ][,
      sigma := sigma
    ]
  )

  ct_trajs <- simulate_cts(params, time_range = time_range, obs_noise = TRUE)

  if (!is.null(sample_density)) {
    ct_trajs <- ct_trajs[,
     .SD[sample %in% sample(.N, sample(sample_density, 1))], by = "id"
    ]
  }

  ct_trajs <- index_by_first_positive(ct_trajs)
  ct_trajs[, onset_time := as.integer(onset_time - t_first_pos)]
  return(ct_trajs)
}
