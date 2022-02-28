  subject_design <- function(formula = ~ 1, data, preds_sd = 0.1, params) {
  subjects <- extract_subjects(data)
  design <- model.matrix(formula, data = subjects)

  out <- list(
    design = design, subjects = subjects, params = params, preds_sd = preds_sd
  )
  return(out)
}


data_to_stan <- function(input_data,
                         ct_model = subject_design(~ 1, input_data),
                         likelihood = TRUE, clod = 40,
                         onsets = TRUE) {


  stan_data <- list(N = input_data[, .N],
                    P = length(unique(input_data$id)),
                    id = input_data[, id],
                    day_rel = input_data[, t],
                    swab_types = length(unique(input_data[, swab_type])) - 1,
                    swab_type = input_data[, swab_type] + 1,
                    ct_value = ifelse(
                      is.na(input_data$ct_value), -99, input_data$ct_value
                    ),
                    pcr_res = input_data[, pcr_res],
                    t_e = 0,
                    c_0 = clod,
                    c_lod = clod,
                    lmean = get_inc_period()$inc_mean_p,
                    lsd = get_inc_period()$inc_sd_p,
                    likelihood = as.numeric(likelihood),
                    preds = ncol(ct_model$design) - 1,
                    preds_sd = ct_model$preds_sd,
                    design = ct_model$design
  )
 if (is.null(input_data$onset_time) | !onsets) {
  stan_data <- c(stan_data, list(
          any_onsets = 0,
          onset_avail = rep(0, stan_data$P),
          onset_time = rep(0, stan_data$P)
        ))
 }else{
  onset_dt <- suppressWarnings(
    input_data[,
    .(onset_time = min(onset_time, na.rm = TRUE)), by = "id"
    ][
      is.infinite(onset_time), onset_time := NA
    ]
  )
  stan_data <- c(stan_data, list(
          any_onsets = 1,
          onset_avail = as.numeric(!is.na(onset_dt$onset_time)),
          onset_time = onset_dt$onset_time %>%
            tidyr::replace_na(0)
        ))
 }

  return(stan_data)
}
