stan_inits <- function(dt) {
  function() {
    inits <- list(
      T_e = purrr::map_dbl(
        1:dt$P,
        ~ truncnorm::rtruncnorm(
          1, a = max(-dt$onset_time[.], 0),
          mean = max(-dt$onset_time[.] + 5, 5), sd = 1
        )
      ),
      c_p_mean = rnorm(1, 0, 1),
      c_p_var = abs(rnorm(1, 0, 0.1)),
      c_p_raw = rnorm(dt$P, 0, 1),
      c_s_mean = rnorm(1, 0, 1),
      c_s_var = abs(rnorm(1, 0, 0.1)),
      c_s_raw = rnorm(dt$P, 0, 1),
      t_p_mean = rnorm(1, 1.61, 0.5),
      t_p_var = abs(rnorm(1, 0, 0.1)),
      t_p_raw = rnorm(dt$P, 0, 1),
      t_s_mean = rnorm(1, 1.61, 0.5),
      t_s_var = abs(rnorm(1, 0, 0.1)),
      t_s_raw = rnorm(dt$P, 0, 1),
      t_lod_mean = rnorm(1, 2.3, 0.5),
      t_lod_var = abs(rnorm(1, 0, 0.1)),
      t_lod_raw = rnorm(dt$P, 0, 1),
      sigma = truncnorm::rtruncnorm(1, a = 0, mean = 3, sd = 0.1)
    )

    if (dt$any_onsets == 1) {
      inits$inc_mean <- rnorm(1, dt$lmean[1], dt$lmean[2])
      inits$inc_sd <- truncnorm::rtruncnorm(
        1, a = 0, mean = dt$lsd[1], sd = dt$lsd[2]
      )
    }
    return(inits)
  }
}