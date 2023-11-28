update_variable_labels <- function(draws, reverse = FALSE) {
  
  draws <- data.table::copy(draws)
  params <- c(
    "c_0", "c_p", "c_s", "t_p", "t_s", "t_lod", "inc_mean", "inc_sd",
    "nat_inc_mean", "nat_inc_sd", "ct_shift", "ct_scale"
  )
  
  clean_params <- c(
    "Ct value at LOD",
    "Ct value at peak",
    "Ct value at switch",
    "Time of peak",
    "Time of switch",
    "Time until PCR-",
    "Incubation period (log mean)",
    "Incubation period (log sd)",
    "Incubation period (mean)",
    "Incubation period (sd)",
    "Ct intercept adjustment",
    "Ct multiplicative adjustment"
  )
  
  if (reverse) {
    params <- rev(params)
    clean_params <- rev(clean_params)
  }
  
  draws <- draws[
    variable %in% params
  ][,
    variable := factor(
      variable,
      levels = params,
      labels = clean_params
    )
  ]
  return(draws[])
}

add_baseline_to_draws <- function(adj_draws, predictor_label, onsets_flag) {
  
  adj_draws <- rbind(
    extract_param_draws(
      draws, onsets = onsets_flag)[, predictor := predictor_label
      ],
    adj_draws
  )[, predictor := factor(predictor)]
  
  return(adj_draws)
  
}

add_regressor_categories <- function(dt) {
  
  dt[predictor %in% c("Omicron (BA.1)", "Omicron (BA.2)", "Delta"),
     regressor_category := "VOC"]
  
  dt[predictor %in% c("Age: 20-34",
                      "Age: 35-49",
                      "Age: 50+"),
     regressor_category := "Age"]
  
  dt[predictor %in% c("3 exposures", "4 exposures",
                      "5+ exposures"), 
     regressor_category := "Number of exposures"]
  
  dt[predictor %in% c("Symptomatic", "Asymptomatic", "Unknown"),
     regressor_category := "Symptom status"]
  
  dt[predictor %like% "Baseline", regressor_category := "baseline"]
  
}

# Add a function mapping labels to each term
update_predictor_labels <- function(dt) {
  dt <- dt[,
           predictor := factor(
             predictor,
             levels =  c(
               "no_vaccines2",
               "symptomsasymptomatic", "symptomsunknown",
               "VOCDelta", "VOCOmicron (BA.2)",
               "time_since_last_dose",
               "no_exposures3", "no_exposures5+",
               "age_group20-34", "age_group50+",
               "swab_typeVTM"
             ),
             labels = c(
               "2 vaccines",
               "Asymptomatic", "Unknown",
               "Delta", "Omicron (BA.2)",
               "Time since last dose",
               "3 exposures", "5+ exposures",
               "Age: 20-34", "Age: 50+",
               "Swab type"
             )
           )
  ]
  return(dt[])
}
