process_data <- function(data_raw) {

  data_proc <- data.table::copy(data_raw)
  
  setnames(
    data_proc, c("ORF1ab", "total infections"), c("ct", "total_infections")
  )
  
  out <- data_proc, swab_date := dmy(swab_date)][
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
    c("N", "number_vaccines (14 days pre ix)"),
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
            "result", "VOC_basis", "centre", "total_infections")
  out[, (facs) := lapply(.SD, factor), .SDcols = facs]
  out[,
    symptoms := factor(
      symptoms, levels = c("0", "1", "unknown"),
      labels = c("asymptomatic", "symptomatic", "unknown")
    )
  ]
  return(out[])
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

  # Drop swab type private
  out <- out[!swab_type %in% "Private"]

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
  return(out)
}

index_by_first_positive <- function(dt) {
  pos_test <- dt[
    pcr_res == 1, .SD[t == min(t)], by = "id"][,
    .(id, t_first_pos = t)
  ]
  dt <- dt[pos_test, on = "id"]
  dt[, t := t - t_first_pos]
  return(dt[])
}
