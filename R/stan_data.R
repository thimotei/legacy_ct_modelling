stan_data_fun <- function(input_data) {
  
  stan_data <- list(N = input_data[, .N], 
                    P = input_data[, uniqueN(ID)],
                    id = input_data[, ID],
                    day_rel = input_data[, t_first_test],
                    ct_value = ifelse(is.na(input_data$ct_adjusted), -99, input_data$ct_scaled),
                    # ct_value = ifelse(is.na(input_data$ct_adjusted), 40, input_data$ct_adjusted),
                    pcr_res = input_data[, as.numeric(result_num)], 
                    # symp_rel = input_data[, max(onset_rel), by = c("ID", "infection_id")][, V1],
                    # te_upper_bound = input_data[, .(te_upper_bound = unique(te_upper_bound)), 
                    #                             by = id] %>% .[, te_upper_bound],
                    t_e = 0,
                    c_0 = (40 - mn)/(mx - mn),
                    c_lod = (40 - mn)/(mx - mn),
                    lmean = get_inc_period()$inc_mean_p[1],
                    lsd = get_inc_period()$inc_sd_p[2]
)
  
  return(stan_data)
}

