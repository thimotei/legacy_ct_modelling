stan_data_simulated <- function(P_arg) {
  
  stan_data <- list(P = P_arg, 
                    t_e = 0,
                    c_0 = (40 - mn)/(mx - mn),
                    c_lod = (40 - mn)/(mx - mn),
                    lmean = EpiNow2::incubation_periods[, mean],
                    lsd = EpiNow2::incubation_periods[, sd])
  
  return(stan_data)
}
