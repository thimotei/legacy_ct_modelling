process_data <- function(dt_raw) {
  
  dt_proc <- data.table::copy(dt_raw)
  
  # Change some column names
  setnames(
    dt_proc, 
    c("ORF1ab", "total infections", "ID", "infection_ID", "number_vaccines"), 
    c("ct_unadjusted", "total_infections", "id", "infection_id", "no_vaccines")
  )
  
  #--- Filtering for only crick tests and tests with ct values + known age
  # Drop private swab types
  out <- dt_proc[centre == "crick" & 
                   ct_unadjusted != "unknown" & 
                   age != "unknown" & 
                   barcode != "LFT" &
                   swab_type != "Private"]
  
  # Drop 1 individual with no vaccine information
  # out[no_vaccines %in% "unknown"]
  out <- out[!(no_vaccines %in% "unknown")]
  
  # Drop 10 individuals with no symptom status (either symptomatic, asymptomatic,
  # or unknown)
  # out[symptoms %in% "unknown"]
  out <- out[!(symptoms %in% "unknown")]
  
  # Drop subject with ID 978 due to large mismatch between onset and positive
  #  tests
  out <- out[id != 978]
  
  # Re-format dates and ct values
  cols <- c("swab_date", "symptom_onset_date",  paste0("date_dose_",1:3))
  out <- out[,(cols) := lapply(.SD, lubridate::dmy), .SDcols = cols
             ][, ct_unadjusted := as.numeric(ct_unadjusted)]
  
  # Fix onset date for this test
  out[barcode == "RLNB620101", symptom_onset_date := lubridate::dmy("18-02-2022")]
  
  # "invalid" and "inconclusive" results are being re-assigned to 
  # positive based on discussions with researchers at the crick
  # that generated the Ct value data
  out[result == "Invalid", result := "Positive"]
  out[result == "Inconclusive", result := "Positive"]
  out[result == "Negative", uncensored := 0]
  out[result == "Positive", uncensored := 1]
  
  # People with 49U barcodes need swab date moving back 1 day
  out <- out[barcode %like% "49U", swab_date := swab_date - 1]
  
  # Set Ct values for negative results
  out <- out[, ct_value := ct_unadjusted
             ][result == "Negative", ct_value := 40]
  
  # Remove entries with NA ct values
  out <- out[!is.na(ct_value)]
  
  #--- counting number of positive swabs (and number of test days) by individual
  no_pos_cts <- out[(result == "Positive" | result == "Inconclusive"), 
                    .(no_pos_results = .N, 
                      ndays = length(unique(swab_date))), 
                    by = c("id", "infection_id")]
  
  # We are interested in number of timepoints, not number of swabs
  # Some people did > 1 swab on one day
  no_pos_cts[, no_pos_results := min(no_pos_results, ndays), 
             c("id", "infection_id")
             ][, ndays := NULL]
  
  out <- merge(
    out, no_pos_cts, by = c("id", "infection_id")
  )
  
  #--- adding time since first positive test by individual
  first_pos_test_date_dt <- out[result == "Positive",
    .(first_pos_test_date = min(swab_date)),
    by = c("id", "infection_id")]
  
  out <- merge(
    out, first_pos_test_date_dt, by = c("id", "infection_id")
  )
  
  # Calculate number of days since first positive test for subsequent tests
  out[,
      time_since_first_pos := as.numeric(
        swab_date - first_pos_test_date, units = "days"
      ),
      by = c("id", "infection_id")
  ]
  
  # naming time since first positive test, also t, just so the same Stan data
  # function works as the simulation study
  out[, t := time_since_first_pos]
  
  # add time at onset
  out[,
      onset_time := as.numeric(symptom_onset_date - first_pos_test_date, units = "days")
  ]
  
  # Deal with people missing symptom status but having an onset date
  out[is.na(symptoms),
      symptoms := ifelse(!is.na(symptom_onset_date), "1", "unknown")]
  
  # Drop infections with gaps between positive tests of more than 60 days
  ids_spurious_gaps <- out[result == "Positive" & 
                             (abs(time_since_first_pos) > 60 | 
                                abs(onset_time) > 60), id]
  
  out <- out[!(id %in% ids_spurious_gaps)]
  
  # Make variables factors
  facs <- c("no_vaccines", "VOC", "swab_type", "dose_1", "dose_2", "dose_3",
            "result", "centre", "total_infections", "symptoms")
  out[, (facs) := lapply(.SD, factor), .SDcols = facs]
  out[,
      symptoms := factor(
        symptoms, levels = c("0", "1", "unknown"),
        labels = c("asymptomatic", "symptomatic", "unknown")
      )
  ]

  # Drop unused factor levels
  out[, (facs) := lapply(.SD, forcats::fct_drop), .SDcols = facs]
  
  # Add age groups
  out[, age.group := cut(x = as.numeric(age), 
                         breaks = c(20, 34, 49, 64, 100), 
                         labels = c("20-34", "35-49", "50-64", "65+"))]
  
  return(out)
}

voc_status_attribution <- function(dt,
                                   group_all = TRUE) {
  
  if(group_all) {
    # assuming that all "variant-like" samples are the same
    dt[VOC %like% "Delta", VOC := "Delta"]
    dt[VOC == "Omicron (BA.1-like)",
       VOC := "Omicron"]
    dt[VOC == "Omicron (BA.2-like)" | VOC == "Omicron-BA2",
       VOC := "BA.2"]
    
    dt[, VOC := factor(VOC,
                       levels = c("Omicron", "Delta", "BA.2"),
                       labels = c("Omicron", "Delta", "BA.2"))]
  }
  else {
    # assuming that all "variant-like" samples are the same
    
    dt[VOC == "Omicron (BA.2-like)" | VOC == "Omicron-BA2",
       VOC := "Omicron (BA.2-like)"]
    dt[VOC == "BA2",
       VOC := "BA.2"]
    
    dt[, VOC := factor(VOC)]
  }
  return(dt[])
}

# Function to postprocess cleaned input data into modelling dataset
subset_data <- function(dt_clean, no_pos_swabs) {
  
  dt_proc <- data.table::copy(dt_clean)
  
  # Filter for number of positive swabs
  out <- dt_proc[no_pos_results >= no_pos_swabs]
  
  # Rename some columns
  out <- out[, t_first_test := time_since_first_pos
    ][, data_id := id
      ][, id := .GRP, by = c("data_id", "infection_id")
        ][, swab_type_num := ifelse(swab_type == "Dry", 0, 1)]
  
  # Set baselines for factors
  out <- out[,
    VOC := forcats::fct_relevel(VOC, "Omicron")
  ][,
    symptoms := forcats::fct_relevel(symptoms, "symptomatic")
  ][,
    no_vaccines := forcats::fct_relevel(no_vaccines, "3")
  ]
  
  # Clean up factor levels
  facs <- c("no_vaccines", "VOC", "swab_type", "dose_1", "dose_2", "dose_3",
            "result", "centre", "total_infections", "symptoms")
  out[, (facs) := lapply(.SD, forcats::fct_drop), .SDcols = facs]
  
  return(out)
}

t_last_dose <- function(dt, imm_delay = 14) {
  
  dt_proc <- data.table::copy(dt)
  
  # Some people have dose dates entered for some tests but not all tests. 
  # We can fill in these missing dose dates so they don't get filtered out 
  dt_proc[, date_dose_1 := ifelse(any(!is.na(date_dose_1)), 
                             max(date_dose_1, na.rm = TRUE),
                             NA), c("id", "infection_id")]
  dt_proc[, date_dose_2 := ifelse(any(!is.na(date_dose_2)), 
                             max(date_dose_2, na.rm = TRUE),
                             NA), c("id", "infection_id")]
  dt_proc[, date_dose_3 := ifelse(any(!is.na(date_dose_3)), 
                             max(date_dose_3, na.rm = TRUE),
                             NA), c("id", "infection_id")]
  
  fz <- function(x, imm_delay){
    m <- unique(x$first_pos_test_date)
    v <- unique(c(x$date_dose_1, x$date_dose_2, x$date_dose_3))
    # remove NA dose dates
    v <- v[!is.na(v)]
    # remove doses after first positive test
    v <- v[v + imm_delay <= m]
    # return time since dose
    out <- as.numeric(m - max(v) - imm_delay)
    return(out)
  }
  
  # rescaling time since last dose data to roughly match the scale of the
  # standard deviation of the effect size parameters
  dt_proc[, time_since_last_dose := fz(.SD, imm_delay = imm_delay) / 400, c("id", "infection_id")]
  
  return(dt_proc)
  
}

index_by_first_positive <- function(dt) {
  pos_test <- dt[
    uncensored == 1, .SD[t == min(t)], by = "id"][,
      .(id, t_first_pos = t)
    ]
  dt <- dt[pos_test, on = "id"]
  dt[, t := t - t_first_pos]
  return(dt[])
}

create_stan_data <- function(dt){
  dt[, new.id := .GRP, c("id", "infection_id")]
  dt[, test.id := 1:.N]
  out <- dt[, .(id = new.id, 
                test.id,
                ct_value,
                test_date = swab_date,
                t,
                onset_date = symptom_onset_date,
                onset_t = onset_time,
                censored = uncensored)]
  return(out)
}
