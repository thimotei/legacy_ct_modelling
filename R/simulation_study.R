prepare_sim_data <- function(P = 50) {
  
  out <- list(
    P = P,
    any_onsets = 1,
    onset_time = rep(0, P),
    c_lod = 40,
    lmean = get_inc_period()$inc_mean_p,
    lsd = get_inc_period()$inc_sd_p,
    preds = 0,
    ct_preds = 0,
    K = 5)
  
  return(out)
}

prepare_ind_params <- function(dt) {
  
  params <- list(
    T_e = purrr::map_dbl(
      1:dt$P,
      ~ truncnorm::rtruncnorm(
        1, a = max(-dt$onset_time[.], 0),
        mean = max(-dt$onset_time[.] + 5, 5), sd = 1
      )
    ),
    c_0 = truncnorm::rtruncnorm(
      1, a = dt$c_lod, mean = dt$c_lod + 10, sd = 5
    ),
    c_p_mean = rnorm(1, 0, 1),
    # c_s_mean = rnorm(1, 0, 1),
    t_p_mean = rnorm(1, 1.61, 0.5),
    # t_s_mean = rnorm(1, 1.61, 0.5),
    t_lod_mean = rnorm(1, 3.22, 0.5),
    ind_var = abs(rnorm(dt$K, 0, 0.2)),
    ind_eta  = matrix(rnorm(dt$P * dt$K, 0, 1), nrow = dt$K, ncol = dt$P),
    sigma = truncnorm::rtruncnorm(1, mean = 0.5, sd = 1, a = 0)
  )
  
  return(params)
}

simulate_ind_params <- function(obs = obs_sim,
                                parameters = prepare_ind_params(obs_sim),
                                time_range = 0:30,
                                sample_density = 2:8,
                                c_p_eff = 1,
                                t_p_eff = 1,
                                t_lod_eff = 1,
                                covariate_description = "baseline") {
  
  inc_mean <- get_inc_period()$inc_mean_p[1]
  inc_sd <- get_inc_period()$inc_sd_p[1]
  
  dt_out <- with(
    parameters,
    data.table::data.table(
      id = 1:obs$P,
      swab_type = "Dry",
      onset_time = rlnorm(obs$P, inc_mean, inc_sd),
      T_e = T_e,
      t_p = exp(t_p_mean*t_p_eff + ind_var[1] * ind_eta[1, ]),
      t_lod = exp(t_lod_mean*t_lod_eff + ind_var[3] * ind_eta[2, ]),
      c_0 = c_0)
    [, c_p := c_0*plogis(c_p_mean*c_p_eff + ind_var[5] * ind_eta[3, ])
    ][, c_lod := obs$c_lod,
    ][, t_lod_abs := t_p + t_lod
    ][, sigma := sigma
    ][, covariate_1 := covariate_description]
    
  )
  
  return(dt_out)
}

sample_from_trajectories <- function(dt_in, 
                                     min_pos = 1,
                                     max_pos = 10,
                                     min_neg = 1,
                                     max_neg = 10) {
  
  dt_pos <- dt_in[ct_value < 40
    ][, .SD[
      sample %in% sample(sample, runif(1, min_pos, max_pos))
      ], by = "id"] 
  
  dt_neg <- dt_in[
    ct_value == 40
  ][, .SD[
    sample %in% sample(sample, runif(1, min_neg, max_neg))
    ], by = "id"] 
  
  dt_out <- rbind(dt_pos, dt_neg)[order(id, sample)]
  
  return(dt_out)
} 