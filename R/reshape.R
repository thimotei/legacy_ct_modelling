reshape_figure_2_panel_a <- function(dt) {
  
  dt_proc <- data.table::copy(dt)
  
  mx_date <- dt_proc[, max(symptom_onset_date, na.rm = TRUE)]
  
  cols_to_keep_vax <- c("id", "infection_id",
                        "date_dose_1",
                        "date_dose_2",
                        "date_dose_3")
  
  dt_vax <- dt_proc[, ..cols_to_keep_vax]
  
  dt_vax[, date_dose_1_end := ymd(date_dose_2 - 1),
         by = c("id", "infection_id")]
  dt_vax[!is.na(date_dose_3), date_dose_2_end := ymd(date_dose_3 - 1),
         by = c("id", "infection_id")]
  dt_vax[is.na(date_dose_3), date_dose_2_end := mx_date,
         by = c("id", "infection_id")]
  dt_vax[, date_dose_3_end := mx_date,
         by = c("id", "infection_id")]
  
  setnames(dt_vax,
           c("date_dose_1",
             "date_dose_2",
             "date_dose_3"),
           c("date_dose_1_start",
             "date_dose_2_start",
             "date_dose_3_start"))
  
  dt_vax_long_start <- melt(dt_vax[, c("id", "infection_id",
                                       "date_dose_1_start",
                                       "date_dose_2_start",
                                       "date_dose_3_start")],
                            measure.vars = c("date_dose_1_start",
                                             "date_dose_2_start",
                                             "date_dose_3_start"),
                            variable.name = "dose",
                            value.name = "start_date") %>%
    unique()
  
  dt_vax_long_end <- melt(dt_vax[, c("id", "infection_id",
                                     "date_dose_1_end",
                                     "date_dose_2_end",
                                     "date_dose_3_end")],
                          measure.vars = c("date_dose_1_end",
                                           "date_dose_2_end",
                                           "date_dose_3_end"),
                          variable.name = "dose",
                          value.name = "end_date") %>%
    unique() 
  
  
  dt_vax_long_start[, dose := str_remove(dose, "date_")][,
    dose := str_remove(dose, "_*_start")
  ]
  dt_vax_long_end[, dose := str_remove(dose, "date_")][,
    dose := str_remove(dose, "_*_end")
  ]
  
  dt_vax_long <- merge.data.table(dt_vax_long_start,
                                  dt_vax_long_end,
                                  by = c("id", "infection_id", "dose"))
  
  
  cols_to_keep_exposures <- c("id", "infection_id",
                              "VOC",
                              "symptom_onset_date",
                              "first_pos_test_date")
  
  dt_exposures <- dt_obs[, ..cols_to_keep_exposures]
  
  dt_exposures_long <- melt(dt_exposures,
                            measure.vars = c("symptom_onset_date",
                                             "first_pos_test_date"),
                            variable.name = "event",
                            value.name = "date") %>% 
    unique()
  
  # adding plot-specific event labels
  dt_exposures_long[event == "first_pos_test_date", event := "Infection"]
  dt_exposures_long[event == "symptom_onset_date", event := "Symptom onset"]
  
  dt_both <- merge.data.table(dt_vax_long,
                              dt_exposures_long,
                              by = c("id", "infection_id"),
                              allow.cartesian = TRUE) %>% 
    melt(., measure.vars = c("VOC", "event"),
         value.name = "VOC/event")
  
  out <- dt_both[, variable := NULL][`VOC/event` != "Infection"][,
  dose := factor(dose, labels = c("Dose 1", "Dose 2", "Dose 3"))][
  order(date), new_id := .GRP, by = c("id")] %>% unique()
  
  return(out)
}
