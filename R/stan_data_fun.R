stan_data_fun <- function(input_data, likelihood = TRUE, clod = 40) {
  
  stan_data <- list(N = input_data[, .N], 
                    P = length(unique(input_data$id)),
                    id = input_data[, id],
                    day_rel = input_data[, t],
                    ct_value = ifelse(
                      is.na(input_data$ct_value), -99, input_data$ct_value
                    ),
                    pcr_res = input_data[, pcr_res],
                    t_e = 0,
                    c_0 = clod,
                    c_lod = clod,
                    lmean = get_inc_period()$inc_mean_p[1],
                    lsd = get_inc_period()$inc_sd_p[2],
                    likelihood = as.numeric(likelihood)
  )
  return(stan_data)
}
