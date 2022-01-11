stan_data_fun <- function(input_data) {
  
  stan_data <- list(N = input_data[, .N], 
       P = input_data[, uniqueN(id)],
       id = input_data[, id],
       day_rel = input_data[, day_rel],
       ct_value = ifelse(is.na(input_data$ct_value), -99, input_data$ct_scaled),
       pcr_res = input_data[, as.numeric(result)], 
       symp_rel = input_data[, .(symp_rel = unique(symp_rel)), by = id][, symp_rel],
       te_upper_bound = input_data[, .(te_upper_bound = unique(te_upper_bound)), 
                                   by = id] %>% .[, te_upper_bound],
       t_e = 0,
       c_0 = (45 - mn)/(mx - mn),
       c_lod = (45 - mn)/(mx - mn),
       lmean = EpiNow2::incubation_periods[, mean],
       lsd = EpiNow2::incubation_periods[, sd])
  
  return(stan_data)
}


