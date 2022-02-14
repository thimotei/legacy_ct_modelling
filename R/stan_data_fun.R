stan_data_fun <- function(input_data) {
  
  stan_data <- list(N = input_data[, .N], 
                    P = length(unique(input_data$id)),
                    id = input_data[, id],
                    day_rel = input_data[, t],
                    ct_value = ifelse(is.na(input_data$ct_value), -99, input_data$ct_value_std),
                    pcr_res = input_data[, pcr_res], 
                    t_e = 0,
                    c_0 = (40 - mn)/(mx - mn),
                    c_lod = (40 - mn)/(mx - mn),
                    lmean = get_inc_period()$inc_mean_p[1],
                    lsd = get_inc_period()$inc_sd_p[2]
  )
  
  return(stan_data)
}
