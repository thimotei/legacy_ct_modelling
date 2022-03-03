process_data <- function(data_raw) {

  setnames(data_raw, "ORF1ab", "ct")
  
  out <- data_raw[, swab_date := dmy(swab_date)][
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
  return(out[])
}

# function to condition full dataset on VOC type and number of positive swabs
# by individual + add unique infection level ID
subset_data <- function(dt_clean_in, voc, no_pos_swabs) {
  dt_out <- dt_clean[VOC %in% voc][,
    t_first_test := as.numeric(swab_date - min(swab_date), units = "days"),
    by = c("id", "infection_id")][
    no_pos_results >= no_pos_swabs][,
    data_id := id][,
    id := .GRP, by = c("data_id", "infection_id")][,
    swab_type_num := as.numeric(!swab_type %in% "Dry")]
  return(dt_out)
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
