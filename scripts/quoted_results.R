source("scripts/setup.R")

# minimum and maximum dates

obs[, min(swab_date, na.rm = TRUE)]
obs[, max(swab_date, na.rm = TRUE)]

# calculating the proportion of symptomatic individuals
n_symptomatic <- obs[, .(symptom_count = uniqueN(id)), "symptoms"]
prop_symp <- n_symptomatic[symptoms == "symptomatic",
                           symptom_count]/n_symptomatic[,sum(symptom_count)]

 
# total number of participants
n_total <- obs[, uniqueN(data_id)]

obs[, uniqueN(data_id), by = VOC]

# numbers of participants stratified by sex and VOC
n_ind_sex_voc <- obs[, .(no = uniqueN(data_id)),
    c("VOC", "sex")][,
  `:=` (VOC = fct_relevel(VOC, "Delta"),
        sex = fct_relevel(sex, "Female")),
  ][order(VOC, sex)]

# just printing table to console
n_ind_sex_voc

# merging totals with subsets for percentage calculation
merge(n_ind_sex_voc,
      n_ind_sex_voc[, .(sum = sum(no)),
                    by = c("VOC")])[, scales::label_percent() (no/sum),
  by = c("VOC", "sex")]

# working out median age and IQR, stratified by VOC
obs[, .(`50%` = as.numeric(median(age)),
        `25%` = as.numeric(median(age)) - IQR(age),
        `75%` = as.numeric(median(age)) + IQR(age)),
    by = "VOC"]

# the number of doses each individual had at the time of their first positive
# test
obs[swab_date <= first_pos_test_date][,
  uniqueN(data_id), c("no_vaccines", "VOC")][,
  VOC := fct_relevel(VOC, "Delta")][order(VOC)]

# the number of doses each individual had at the time of their first positive
# test minus 14 days (the amount of time required for immunity to build up)

# for 2 doses
obs[first_pos_test_date < date_dose_2 + 14][, 
  uniqueN(data_id), by = "VOC"]

# for 3 doses
obs[first_pos_test_date < date_dose_3 + 14][, 
  uniqueN(data_id), by = "VOC"]

# numbers of individuals stratified by type of dose
# dose 1
obs[, uniqueN(data_id), c("dose_1", "VOC")]

# dose 2
obs[, uniqueN(data_id), c("dose_2", "VOC")]

# dose 3
obs[, uniqueN(data_id), c("dose_3", "VOC")]

# number of individuals stratified by total number of exposures
# (at time of first positive test)
obs[, uniqueN(data_id), c("no_exposures", "VOC")]
 
# median (and IQR) of the time since the last dose prior to infection
obs[, .(as.numeric(median(t_since_last_dose)),
        as.numeric(median(t_since_last_dose)) - IQR(t_since_last_dose),
        as.numeric(median(t_since_last_dose)) + IQR(t_since_last_dose)),
    by = VOC]

obs[, .(quantile(t_since_last_dose, 0.5),
        quantile(t_since_last_dose, 0.25),
        quantile(t_since_last_dose, 0.75)),
    by = VOC]


.obs[, .(IQR(t_since_last_dose)),
    by = VOC]


# either run inference script or load in saved fit
# fit <- readRDS("outputs/fits/fit_full.rds")

cols <- c("variable", "predictor", "median", "lo95", "hi95")
cols_fun <- c("median", "lo95", "hi95")

effect_size_summary_natural[regressor_category == "VOC"][order(variable, predictor)][, ..cols][ , (cols_fun) := lapply(.SD, signif, 2), .SDcols = cols_fun][]
effect_size_summary_natural[regressor_category == "Symptom status"][order(variable, predictor)][, ..cols][ , (cols_fun) := lapply(.SD, signif, 2), .SDcols = cols_fun][]
effect_size_summary_natural[regressor_category == "Number of exposures"][order(variable, predictor)][, ..cols][ , (cols_fun) := lapply(.SD, signif, 2), .SDcols = cols_fun][]
effect_size_summary_natural[regressor_category == "Age"][order(variable, predictor)][, ..cols][ , (cols_fun) := lapply(.SD, signif, 2), .SDcols = cols_fun][]

# recalculating the population-level Ct curves summarised, without censoring of Ct values at LOD
pop_ct_draws_sum <- ct_pp_sum_data(pop_ct_draws_cens)

summarise_positivity_times(pop_ct_draws_sum, ct_threshold = 37)
summarise_positivity_times(pop_ct_draws_sum, ct_threshold = 30)

fig_4_sum <- dt_figure_4[, .(t_me = quantile(t, 0.5),
                t_lo95  = quantile(t, 0.0025),
                t_hi95 = quantile(t, 0.975),
                inc_me = signif(quantile(inc_mean_nat, 0.5), 2),
                inc_lo95  = signif(quantile(inc_mean_nat, 0.0025), 2),
                inc_hi95 = signif(quantile(inc_mean_nat, 0.975), 2)),
            by = c("regressor_category", "predictor")][order(regressor_category)]

fig_4_sum[regressor_category == "VOC"]
fig_4_sum[regressor_category == "Age"]




