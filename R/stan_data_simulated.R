stan_data_simulated <- function(P_arg, t_max_arg) {
  
  stan_data <- list(P = P_arg,
                    t_max = t_max_arg,
                    t_e = 0,
                    c_0 = (40 - mn)/(mx - mn),
                    c_lod = (40 - mn)/(mx - mn),
                    lmean = get_inc_period()$inc_mean_p[1],
                    lsd = get_inc_period()$inc_sd_p[2]
)
  return(stan_data)
}
