# calculate_ct_threshold <- function(dt,
#                                    ct_threshold,
#                                    trim_flag) {
#   
#   dt_proc <- data.table::copy(dt)
#   # out <- dt_proc[, .SD[which.min(abs(ct_value - ct_threshold)) & t < t_p],
#   #                by = c(".draw", "id", "predictor", "regressor_category")]
#   # 
#   out <- dt_proc[, .SD[t < t_p],
#     by = c(".draw", "id", "predictor", "regressor_category")][, 
#     .SD[which.min(abs(ct_value - ct_threshold))],
#     by = c(".draw", "id", "predictor", "regressor_category")]
#   
#   if(trim_flag == TRUE) {
#     out <- out[, c("t", "ct_value", 
#                    ".draw", "id", "predictor",
#                    "regressor_category")]
#   }
#   
#   return(out)
# }

calculate_t_ct_threshold <- function(dt_in, 
                                     ct_threshold) {
  
  pre_peak <- dt_in[t < t_p, .SD[which.min(abs(ct_value - ct_threshold))],
                              by = c(".draw",
                                     "inc_mean_nat",
                                     "predictor",
                                     "regressor_category")
  ][, c(".draw", 
        "predictor", 
        "regressor_category", 
        "t")
  ]
  
  pre_peak[, direction := "up"]
  
  post_peak <- dt_in[t > t_p, .SD[which.min(abs(ct_value - ct_threshold))],
                               by = c(".draw",
                                      "inc_mean_nat",
                                      "predictor",
                                      "regressor_category")
  ][, c(".draw", 
        "predictor",
        "regressor_category", 
        "t")
  ]
  
  post_peak[, direction := "down"]
  
  dt_out <- rbind(pre_peak, post_peak)
  
  return(dt_out)
}
