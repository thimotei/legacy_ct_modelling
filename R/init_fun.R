init_fun <- function(chain_id) {
  
  out <- list(
    T_e = rnorm(stan_data_simulated$P, 0 , 1),
    
    t_p_mean = runif(1, 0, 1), # rnorm(n = 1, mean = log(5), sd = 0.1)
    t_p_var = 1, # rtruncnorm(n = 1, mean = 0, sd = 0.25, a = 0)
    t_p_raw = rep(0, stan_data_simulated$P), # rnorm(n = P, mean = 0, sd = 1)
    
    t_s_mean = runif(1, 0, 1), # rnorm(n = 1, mean = log(5), sd = 0.1)
    t_s_var = 1, # rtruncnorm(n = 1, mean = 0, sd = 0.25, a = 0)
    t_s_raw = rep(0, stan_data_simulated$P), # rnorm(n = P, mean = 0, sd = 1)
    
    t_lod_mean = runif(1, 0, 1), # rnorm(n = 1, mean = log(10), sd = 0.1)
    t_lod_var = 1, # rtruncnorm(n = 1, mean = 0, sd = 0.25, a = 0)
    t_lod_raw = rep(0, stan_data_simulated$P), # rnorm(n = P, mean = 0, sd = 1)
    
    c_p_mean = 0, # rnorm(1, 0, 0.25)
    c_p_var = 1, # rtruncnorm(n = 1, mean = 1, sd = 0.25, a = 0)
    c_p_raw = rep(0, stan_data_simulated$P), # rnorm(P, 0 , 0.1)
    
    c_s_mean = 0, # rnorm(1, 0, 1)
    c_s_var = 1, # rtruncnorm(n = 1, mean = 0, sd = 1, a = 0)
    c_s_raw = rep(0, stan_data_simulated$P), # rnorm(P, 0 , 0.1)
    
    sigma_obs = 1
  )
  
  return(out)
}