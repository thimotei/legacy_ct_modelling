process_data <- function(dt,
                         struct_arg = "long") {
  
  dt <- fread("data/raw_data.csv")
  
  dt_proc <- data.table::copy(dt)
  
  # Change some column names
  setnames(
    dt_proc, 
    c("ORF1ab", "total infections", "ID", "infection_ID",
      "number_vaccines", "N gene Ct", "S gene Ct"), 
    c("ct_unadjusted", "total_infections", "id", "infection_id",
      "no_vaccines", "ct_n_gene", "ct_s_gene")
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
  # out <- out[!(symptoms %in% "unknown")]
  
  # Drop subject with ID 978 due to large mismatch between onset and positive
  # tests
  out <- out[id != 978]
  
  # Re-format dates and ct values
  cols <- c("swab_date", "symptom_onset_date", paste0("date_dose_", 1:4))
  out <- out[,
             (cols) := lapply(.SD, lubridate::dmy), .SDcols = cols][,
                                                                    ct_unadjusted := as.numeric(ct_unadjusted)][, 
                                                                                                                ct_n_gene := as.numeric(ct_n_gene)][,
                                                                                                                                                    ct_s_gene := as.numeric(ct_s_gene)]
  
  # Fix onset date for this infection episode
  out[barcode == "RLNB620101",
      symptom_onset_date := lubridate::dmy("18-02-2022")]
  
  # Fix onset date for this infection episode
  out[id == 857,
      symptom_onset_date := lubridate::dmy("19-02-2022")]
  
  # fix symptom status for the same individual
  out[id == 857,
      symptoms := "1"]
  
  # fixing the number of doses for an individual with some incorrect data
  out[id == 135, no_vaccines := 3]
  out[id == 135, dose_3 := "BNT162b2"]
  
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
  out <- out[!(is.na(ct_value) & is.na(ct_n_gene) & is.na(ct_s_gene))]
  
  # Determining positivity or negativity using Ct values at 
  # any one of the three gene targets
  out <- out[ct_value < 37 | ct_n_gene < 37 | ct_s_gene < 37,
             result := "Positive"]
  out <- out[ct_value >= 37 | ct_n_gene >= 37 | ct_s_gene >= 37,
             result := "Negative"]
  
  # Counting number of positive swabs (and number of test days) by individual
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
      t_since_first_pos := as.numeric(
        swab_date - first_pos_test_date, units = "days"
      ),
      by = c("id", "infection_id")
  ]
  
  # naming time since first positive test, also t, just so the same Stan data
  # function works as the simulation study
  out[, t := t_since_first_pos]
  
  # add time at onset
  out[,
      onset_time := as.numeric(symptom_onset_date - first_pos_test_date,
                               units = "days")
  ]
  
  # calculating number of exposures as the sum of number of vaccines and
  # number of infections by individual
  out[, no_exposures := as.numeric(no_vaccines) + as.numeric(total_infections),
      by = id]
  
  # Add symptomatic status to individuals with missing status but onset date
  # as per advice from Crick collaborators
  out[is.na(symptoms),
      symptoms := ifelse(!is.na(symptom_onset_date), "1", "unknown")]
  
  # Drop infections with gaps between positive tests of more than 60 days
  ids_spurious_gaps <- out[result == "Positive" & 
                             (abs(t_since_first_pos) > 60 | 
                                abs(onset_time) > 60), id]
  out <- out[!(id %in% ids_spurious_gaps)]
  
  # Make variables factors
  facs <- c("no_vaccines", "VOC", "swab_type", "dose_1", "dose_2", "dose_3",
            "result", "centre", "total_infections", "symptoms", "no_exposures")
  out[, (facs) := lapply(.SD, factor), .SDcols = facs]
  
  # combining 5 and greater than 5 exposures,
  # as only a few individuals have more than 5
  out[no_exposures == "5" | no_exposures == "6",
      no_exposures := factor("5+"),
      by = id][,
               no_exposures := forcats::fct_drop(no_exposures)]
  
  # re-assign those with unknown symptom status to symptomatic (for the time-being,
  # based on a discussion with Crick collaborators)
  out[,
      symptoms := factor(symptoms,
                         levels = c("0", "1", "unknown"),
                         labels = c("asymptomatic", "symptomatic", "symptomatic")
      )
  ]
  
  # Drop unused factor levels
  out[, (facs) := lapply(.SD, forcats::fct_drop), .SDcols = facs]
  
  # Add age groups
  out[, age_group := cut(x = as.numeric(age), 
                         breaks = c(20, 34, 49, 100), 
                         labels = c("20-34", "35-49", "50+"))]
  
  
  if(struct_arg == "long") {
    # melting the Ct values for the three different gene targets, so we can
    # use all three in the inference
    out <- melt(out,
                measure.vars = c("ct_value", "ct_n_gene", "ct_s_gene"),
                variable.name = "ct_type", 
                value.name = "ct_value")[!is.na(ct_value)]
  }
  
  return(out)
}

voc_status_attribution <- function(dt) {
  
  out <- data.table::copy(dt)
  
  # assuming that all "variant-like" samples are the same
  out[VOC %like% "Delta", VOC := "Delta"]
  out[VOC == "Omicron (BA.1-like)",
      VOC := "Omicron (BA.1)"]
  out[VOC == "Omicron (BA.2-like)" | VOC == "Omicron-BA2",
      VOC := "Omicron (BA.2)"]
  
  out[, VOC := factor(VOC,
                      levels = c("Omicron (BA.1)",
                                 "Delta",
                                 "Omicron (BA.2)"),
                      labels = c("Omicron (BA.1)",
                                 "Delta",
                                 "Omicron (BA.2)"))]
  
  
  # removing any data without a VOC attribution
  out <- out[!is.na(VOC)]
  
  return(out[])
}

# Function to postprocess cleaned input data into modelling dataset
subset_data <- function(dt, no_pos_swabs, first_pos_min) {
  
  dt_proc <- data.table::copy(dt)
  
  # Filter for number of positive swabs
  out <- dt_proc[no_pos_results >= no_pos_swabs]
  
  # Rename some columns
  out <- out[, t_first_test := t_since_first_pos
  ][, data_id := .GRP, by = id
  ][, id := .GRP, by = c("data_id", "infection_id")
  ][, swab_type_num := ifelse(swab_type == "Dry", 0, 1)]
  
  # Set baselines for factors
  out <- out[,
             VOC := forcats::fct_relevel(VOC, "Omicron (BA.1)")
  ][,
    symptoms := forcats::fct_relevel(symptoms, "symptomatic")
  ][,
    no_exposures := forcats::fct_relevel(no_exposures, "4")
  ][,
    no_vaccines := forcats::fct_relevel(no_vaccines, "3")
  ][,
    age_group := forcats::fct_relevel(age_group, "35-49")]
  
  # Clean up factor levels
  facs <- c("no_vaccines", "VOC", "swab_type", "dose_1", "dose_2", "dose_3",
            "result", "centre", "total_infections", "symptoms", "no_exposures")
  
  out[, (facs) := lapply(.SD, forcats::fct_drop), .SDcols = facs]
  
  out <- out[!t_since_first_pos < first_pos_min]
  
  return(out)
}

t_first_inf <- function(dt, imm_delay = 14) {
  
  dt_proc <- data.table::copy(dt)
  
  dt_proc <- dt_proc[, t_first_inf_by_id := min(first_pos_test_date),
                     by = "id"
  ][, t_first_inf_to_cur_inf := as.numeric(first_pos_test_date - t_first_inf_by_id),
    by = c("id", "infection_id")]
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
  dt_proc[, date_dose_4 := ifelse(any(!is.na(date_dose_4)),
                                  max(date_dose_4, na.rm = TRUE),
                                  NA), c("id", "infection_id")]
  
  fz <- function(x, imm_delay){
    m <- unique(x$first_pos_test_date)
    v <- unique(c(x$date_dose_1, x$date_dose_2,
                  x$date_dose_3, x$date_dose_4))
    # remove NA dose dates
    v <- v[!is.na(v)]
    # remove doses after first positive test
    v <- v[v + imm_delay <= m]
    # return time since dose
    out <- as.numeric(m - max(v) - imm_delay)
    return(out)
  }
  
  # calculating time since last dose in days
  dt_proc[, t_since_last_dose := fz(.SD, imm_delay = imm_delay),
          c("id", "infection_id")]
  
  return(dt_proc)
  
}

t_last_exposure <- function(dt, imm_delay = 14) {
  
  dt_proc <- data.table::copy(dt)
  
  dt_proc <- t_last_dose(dt_proc, imm_delay)
  dt_proc <- t_first_inf(dt_proc, imm_delay)
  
  dt_proc[t_first_inf_to_cur_inf != 0,
          t_since_last_exposure := min(pmin(t_first_inf_to_cur_inf,
                                            t_since_last_dose)),
          by = c("id", "infection_id")]
  
  dt_proc[t_first_inf_to_cur_inf == 0,
          t_since_last_exposure := t_since_last_dose, 
          by = c("id", "infection_id")]
  
  # rescaling using min-max feature rescaling, normalising
  # the time since last event units to a [0, 1] range
  # used here for rescaling the time since the last exposure
  mn <- dt_proc[, min(t_since_last_exposure)]
  mx <- dt_proc[, max(t_since_last_exposure)]
  
  dt_proc[,
          t_since_last_exposure := (t_since_last_exposure - mn)/(mx - mn),
          c("id", "infection_id")
  ]
  
  return(dt_proc)
}

index_by_first_positive <- function(dt) {
  
  out <- data.table::copy(dt)
  
  pos_test <- out[
    uncensored == 1, .SD[t == min(t)], by = "id"][,
                                                  .(id, t_first_pos = t)
    ]
  
  out <- out[pos_test, on = "id"]
  out[, t := t - t_first_pos]
  
  return(out[])
}

create_stan_data <- function(dt){
  
  out <- data.table::copy(dt)
  
  out[, new.id := .GRP, c("id", "infection_id")]
  out[, test.id := 1:.N]
  
  out <- out[, .(id = new.id, 
                 test.id,
                 ct_value,
                 test_date = swab_date,
                 t,
                 onset_date = symptom_onset_date,
                 onset_t = onset_time,
                 censored = uncensored)]
  
  return(out)
}