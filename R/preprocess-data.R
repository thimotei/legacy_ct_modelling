process_data <- function(dt_raw) {
  
  dt_proc <- data.table::copy(dt_raw)
  
  setnames(
    dt_proc, c("ORF1ab", "total infections"), c("ct", "total_infections")
  )
  
  out <- dt_proc[, swab_date := dmy(swab_date)][
    barcode %like% "49U", swab_date := swab_date - 1][
    symptom_onset_date == "unknown", symptom_onset_date := NA][,
    symptom_onset_date := dmy(symptom_onset_date)][,
    date_dose_1 := dmy(date_dose_1)][,
    date_dose_2 := dmy(date_dose_2)][,
    date_dose_3 := dmy(date_dose_3)][
    ct != "unknown"][,
    ct := as.numeric(ct)][
    result == "Negative", ct := 40][
    swab_type == "VTM" & result == "Negative", ct_adjusted := 40][,
    ct_adjusted := ct]
  
  #--- conditioning on only tests performed in the Crick, leaving out
  #--- UCLH swabs for now
  out <- out[centre == "crick"] # nolint
  
  #--- counting number of positive swabs by individual
  no_pos_cts <- out[
    (result == "Positive" | result == "Inconclusive") &
      is.na(ct_adjusted) == FALSE][,
                                   .N, by = c("ID", "infection_ID")]
  
  out <- merge(
    out, no_pos_cts, by = c("ID", "infection_ID")
  )
  
  setnames(
    out,
    c("N", "number_vaccines"),
    c("no_pos_results", "no_vaccines")
  )
  out[result == "Inconclusive"]
  
  setnames(
    out,
    c("ID", "infection_ID", "ct", "ct_adjusted"),
    c("id", "infection_id", "ct_unadjusted", "ct_value")
  )
  
  #--- conditioning out observations with NA Ct values (how can we use these?)
  out <- out[!is.na(ct_unadjusted)]
  
  #--- adding time since first positive test by individual
  first_pos_test_date_dt <- out[
    result == "Positive",
    .(first_pos_test_date = min(swab_date)),
    by = c("id", "infection_id")
  ]
  
  out <- merge(
    out, first_pos_test_date_dt, by = c("id", "infection_id")
  )
  
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
      onset_time := as.numeric(symptom_onset_date - first_pos_test_date)
  ]
  
  # Deal with people missing symptom status but having an onset date
  out[!symptoms %in% c("0", "1", "unknown"),
      symptoms := ifelse(!is.na(symptom_onset_date), "1", "unknown")
  ]
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
  
  # "invalid" and "inconclusive" both have many
  # observations consistent with positive Ct values, keeping them for now
  # will discuss with Crick partners about this
  out[result == "Invalid", result := "Positive"]
  out[result == "Inconclusive", result := "Positive"]
  out[result == "Negative", uncensored := 0]
  out[result == "Positive", uncensored := 1]
  
  # Drop infections with gaps between positive tests of more than 60 days
  ids_spurious_gaps <- out[,
    .(spurious = any(
      abs(time_since_first_pos) > 60) || any(abs(onset_time) > 60)),
    by = c("id", "infection_id")]
  ids_spurious_gaps[spurious == TRUE][]
  out <- out[ids_spurious_gaps, on = "id"]
  out <- out[spurious == FALSE | is.na(spurious)]
  
  # Drop subject with ID 978 due to large mismatch between onset and positive
  #  tests
  out <- out[id != 978]
  
  # Drop 1 individual with no vaccine information
  out[no_vaccines %in% "unknown"]
  out <- out[!no_vaccines %in% "unknown"]
  
  # Drop 2 individuals with no symptom status (either symptomatic, asymptomatic,
  # or unknown)
  out[symptoms %in% "unknown"]
  out <- out[!symptoms %in% "unknown"]
  # Drop swab type private
  out <- out[!swab_type %in% "Private"]
  out[, (facs) := lapply(.SD, forcats::fct_drop), .SDcols = facs]
  return(out[])
}

voc_status_attribution <- function(dt,
                                   group_all = TRUE) {
  
  if(group_all) {
    # assuming that all "variant-like" samples are the same
    dt[VOC == "Delta (B.1.617.2-like)" | VOC == "Delta (AY.4-like)",
       VOC := "Delta"]
    dt[VOC == "Omicron (BA.1-like)",
       VOC := "Omicron"]
    dt[VOC == "Omicron (BA.2-like)" | VOC == "Omicron-BA2",
       VOC := "BA.2"]
    dt[VOC == "BA2",
       VOC := "BA.2"]
    
    dt[, VOC := factor(VOC,
                       levels = c("Omicron", "Delta", "BA.2", "unknown"),
                       labels = c("Omicron", "Delta", "BA.2", "Unknown"))]
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
  
  out <- dt_proc[,
    t_first_test := as.numeric(swab_date - min(swab_date), units = "days"),
    by = c("id", "infection_id")][
    no_pos_results >= no_pos_swabs][,
    data_id := id][,
    id := .GRP, by = c("data_id", "infection_id")][,
    swab_type_num := as.numeric(!swab_type %in% "Dry")]
  
  # Assume potential BA.2 are BA.2
  out <- out[VOC %in% "?BA2", VOC := "BA2"][,
    VOC := forcats::fct_drop(VOC)
  ]
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
  return(out[])
}

time_since_last_dose <- function(dt, imm_delay = 14) {
  
  dt[!is.na(dose_1) & !is.na(date_dose_1) & 
       date_dose_2 + imm_delay > first_pos_test_date,
       time_since_last_dose := as.numeric(first_pos_test_date - date_dose_1, units = "days"),
     by = c("id", "infection_id")]
  
  dt[!is.na(dose_2) & !is.na(date_dose_2) & 
       date_dose_2 + imm_delay <= first_pos_test_date,
     time_since_last_dose := as.numeric(first_pos_test_date - date_dose_2, units = "days"),
     by = c("id", "infection_id")]
  
  dt[!is.na(dose_3) & !is.na(date_dose_3) &
       date_dose_3 + imm_delay <= first_pos_test_date, 
     time_since_last_dose := as.numeric(first_pos_test_date - date_dose_3, units = "days"),
     by = c("id", "infection_id")]
  
  dt <- dt[!(is.na(date_dose_1) | is.na(date_dose_2) | is.na(date_dose_3))]
  
  # rescaling time since last dose data to roughly match the scale of the
  # standard deviation of the effect size parameters
  dt[, time_since_last_dose := time_since_last_dose/400]
  
  return(dt)
  
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
