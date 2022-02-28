data_to_stan <- function(input_data, form = ~ + 1,
                         coeff_sd = 1,
                         likelihood = TRUE, clod = 40,
                         onsets = TRUE) {

  design <- model.matrix(form, data = input_data)

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
                    preds = ncol(design) - 1,
                    preds_sd = coeff_sd,
                    X = design
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
