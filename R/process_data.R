process_data <- function(data_raw, adj_params) {
    
  #--- dry vs wet swab adjustment. We adjust the VTM swabs downwards
  adjustment <- function(alpha, beta, x) {
    alpha + beta * x
  }

  data_out <- data_raw[, swab_date := dmy(swab_date)] %>% 
    .[barcode %like% "49U", swab_date := swab_date - 1] %>% 
    .[symptom_onset_date == "unknown", symptom_onset_date := NA] %>% 
    .[, symptom_onset_date := dmy(symptom_onset_date)] %>% 
    .[, date_dose_1 := dmy(date_dose_1)] %>% 
    .[, date_dose_2 := dmy(date_dose_2)] %>% 
    .[, date_dose_3 := dmy(date_dose_3)] %>% 
    .[ct != "unknown"] %>% 
    .[, ct := as.numeric(ct)] %>% 
    .[result == "Negative", ct := 40] %>% 
    .[swab_type == "VTM" & result == "Negative", ct_adjusted := 40] %>%
    .[swab_type == "VTM" & result != "Negative",
      ct_adjusted := adjustment(
        adj_params[param == "alpha", me], adj_params[param == "beta", me], ct
      )
    ] %>%
    .[swab_type != "VTM", ct_adjusted := ct]
  
  #--- removing duplicate VTM vs Dry swabs
  barcodes.to.remove <-  data_out[swab_type == "Dry" | swab_type == "VTM"] %>% 
    .[, total_obs_by_date := .N, by = c("ID", "swab_date")] %>% 
    .[total_obs_by_date > 1, c("ID", "swab_date", "barcode", "swab_type")] %>% 
    .[order(ID, swab_date)] %>% 
    .[swab_type == "Dry", barcode]
  
  data_out_no_dups <- data_out[!barcode %in% barcodes.to.remove] 
  
  #--- conditioning on only tests performed in the Crick, leaving out
  #--- UCLH swabs for now
  data_out_no_dups <- data_out_no_dups[centre == "crick"]
  
  #--- counting number of positive swabs by individual
  no_pos_cts <- data_out_no_dups[(result == "Positive" |
                                  result == "Inconclusive") & 
                                  is.na(ct_adjusted) == FALSE] %>% 
    .[, .N, by = c("ID", "infection_ID")]
  
  data_out_no_dups <- merge.data.table(data_out_no_dups,
                                       no_pos_cts, 
                                       by = c("ID", "infection_ID")) %>%  
    setnames(., 
             c("N", "number_vaccines (14 days pre ix)"), 
             c("no_pos_results", "no_vaccines"))
  
  data_out_no_dups[result == "Inconclusive"]
  
  setnames(data_out_no_dups,
           c("ID", "infection_ID", "ct", "ct_adjusted"),
           c("id", "infection_id", "ct_unadjusted", "ct_value"))
  
  #--- conditioning out observations with NA Ct values (how can we use these?)
  data_out_no_dups <- data_out_no_dups[is.na(ct_unadjusted) == FALSE]
  
  #--- adding time since first positive test by individual
  first_pos_test_date_dt <- data_out_no_dups[result == "Positive",
                                     .(first_pos_test_date = min(swab_date)),
                                     by = c("id", "infection_id")]

  data_out_no_dups <- merge.data.table(data_out_no_dups,
                   first_pos_test_date_dt,
                   by = c("id", "infection_id"))

  data_out_no_dups[, time_since_first_pos := as.numeric(swab_date - first_pos_test_date, units = "days"),
                   by = c("id", "infection_id")]
  
  # naming time since first positive test, also t, just so the same Stan data
  # function works as the simulation study
  data_out_no_dups[, t := time_since_first_pos]
    
  return(data_out_no_dups)
}
