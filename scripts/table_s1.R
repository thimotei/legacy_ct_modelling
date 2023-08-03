# isolating all individuals with no infection episodes before joining the study
# remember this does not necessarily mean they have not been infected. this
# just means they hadn't had a detected infection before the start of the study
load("~/lshtm/legacy_correlates/data/Legacy_DataAnnotatedDateSeries_2023-05-04.RData")

# checking everyone in the Ct study is in the full dataset 
ct_study_ids <- obs[, unique(data_id)]

dt_legacy_raw[as.numeric(elig_study_id) %in% ct_study_ids, uniqueN(elig_study_id)]
# three IDs are missing from the full legacy study dataset


dt_legacy_raw <- data.table(chrono.df)

dt_legacy_no_inf_before_start <- dt_legacy_raw[
  episode_number == 1
  ][, .SD[episode_start > enrolment_date], 
    by = elig_study_id][, .(id = as.numeric(elig_study_id), episode_start, enrolment_date)
  ] |> unique()

dt_legacy_inf_before_start <- dt_legacy_raw[episode_number == 1][,
  .SD[episode_start <= enrolment_date], 
  by = elig_study_id][, .(id = as.numeric(elig_study_id), episode_start, enrolment_date)
  ] |> unique()

dt_legacy_raw[episode_number == 1, uniqueN(elig_study_id)]

ids_no_inf_before_start <- dt_legacy_no_inf_before_start[, id]
ids_inf_before_start <- dt_legacy_inf_before_start[, id]

dt_ct_no_inf_before_start <- obs[data_id %in% ids_no_inf_before_start]
dt_ct_inf_before_start <- obs[data_id %in% ids_inf_before_start]

p1 <- dt_ct_no_inf_before_start[result == "Positive"][, VOC := fct_relevel(VOC, "Delta")] |> 
  ggplot(aes(x = VOC, y = ct_value, colour = VOC)) + 
  # geom_point() +
  geom_violin(alpha = 0.1) + 
  geom_jitter(height = 0, width = 0.1) +
  theme_minimal() + 
  scale_y_reverse(limits = c(40, 10)) 
  
p2 <- dt_ct_inf_before_start[result == "Positive"][, VOC := fct_relevel(VOC, "Delta")] |> 
  ggplot(aes(x = VOC, y = ct_value, colour = VOC)) + 
  # geom_point() +
  geom_violin(alpha = 0.1) + 
  geom_jitter(height = 0, width = 0.1) +
  theme_minimal() + 
  scale_y_reverse(limits = c(40, 10)) 

p1/p2 + plot_layout(guides = "collect")

dt_ct_no_inf_before_start[, uniqueN(data_id)]
dt_ct_inf_before_start[, uniqueN(data_id)]
obs[, uniqueN(data_id)]

# need to make a new table with the following:
# number of points available per patient
# delay between first viral load data and symptom onset
# observed peak viral load
# time to clearance

#--- Summary statistics for infection naive individuals
# number of points per patient
dt_ct_no_inf_before_start[, .N, by = c("data_id", "VOC", "result")
][, .(me = signif(quantile(N, 0.5), 2), 
      lo = signif(quantile(N, 0.5) - IQR(N)/2, 2),
      hi = signif(quantile(N, 0.5) + IQR(N)/2, 2)), by = VOC]

dt_ct_no_inf_before_start[, .N, by = c("data_id", "VOC", "result")
][, median(N), by = VOC]

# delay between first viral load data and symptom onset
dt_onset_times_no_inf <- dt_ct_no_inf_before_start[, .(onset_time), by = c("data_id", "VOC")
                                                   ][!is.na(onset_time)] |>
  unique()

dt_onset_times_no_inf[, .(
  me = signif(quantile(onset_time, 0.5), 2), 
  lo = signif(quantile(onset_time, 0.5) - IQR(onset_time)/2, 2),
  hi = signif(quantile(onset_time, 0.5) + IQR(onset_time)/2, 2)), by = VOC]


# observed peak viral load
dt_ct_no_inf_before_start[, .SD[which.min(ct_value)], by = id
][, .(me = signif(median(ct_value), 3), 
      lo = signif(median(ct_value) - IQR(ct_value)/2, 3),
      hi = signif(median(ct_value) + IQR(ct_value)/2, 3)), by = VOC]


# time to clearance
dt_ct_no_inf_before_start[result == "Positive", .SD[which.max(t)], by = id
][, .(me = signif(median(t), 2), 
      lo = signif(median(t) - IQR(t)/2, 2),
      hi = signif(median(t) + IQR(t)/2, 2)), by = VOC]

#--- Summary statistics for previously infected individuals
# number of points per patient
dt_ct_inf_before_start[, .N, by = c("data_id", "VOC", "result")
][, .(me = signif(quantile(N, 0.5), 2), 
      lo = signif(quantile(N, 0.5) - IQR(N)/2, 2),
      hi = signif(quantile(N, 0.5) + IQR(N)/2, 2)), by = VOC]

# delay between first viral load data and symptom onset
dt_onset_times_inf <- dt_ct_inf_before_start[, .(onset_time), by = c("data_id", "VOC")
][!is.na(onset_time)] |>
  unique()

dt_onset_times_inf[, .(
  me = signif(quantile(onset_time, 0.5), 2), 
  lo = signif(quantile(onset_time, 0.5) - IQR(onset_time)/2, 2),
  hi = signif(quantile(onset_time, 0.5) + IQR(onset_time)/2, 2)), by = VOC]

# observed peak viral load
dt_ct_inf_before_start[, .SD[which.min(ct_value)], by = id
][, .(me = signif(median(ct_value), 3), 
      lo = signif(median(ct_value) - IQR(ct_value)/2, 3),
      hi = signif(median(ct_value) + IQR(ct_value)/2, 3)), by = VOC]
# time to clearance
dt_ct_inf_before_start[result == "Positive", .SD[which.max(t)], by = id
][, .(me = signif(median(t), 2), 
      lo = signif(median(t) - IQR(t)/2, 2),
      hi = signif(median(t) + IQR(t)/2, 2)), by = VOC]

#--- Summary statistics for whole cohort
# number of points per patient
obs[, .N, by = c("data_id", "VOC", "result")
][, .(me = signif(quantile(N, 0.5), 2), 
      lo = signif(quantile(N, 0.5) - IQR(N)/2, 2),
      hi = signif(quantile(N, 0.5) + IQR(N)/2, 2)), by = VOC]

# delay between first viral load data and symptom onset
dt_onset_times_all <- obs[, .(onset_time), by = c("data_id", "VOC")
][!is.na(onset_time)] |>
  unique()

dt_onset_times_all[, .(
  me = signif(quantile(onset_time, 0.5), 2), 
  lo = signif(quantile(onset_time, 0.5) - IQR(onset_time)/2, 2),
  hi = signif(quantile(onset_time, 0.5) + IQR(onset_time)/2, 2)), by = VOC]

# observed peak viral load
obs[, .SD[which.min(ct_value)], by = id
][, .(me = signif(median(ct_value), 3), 
      lo = signif(median(ct_value) - IQR(ct_value)/2, 3),
      hi = signif(median(ct_value) + IQR(ct_value)/2, 3)), by = VOC]

# time to clearance
obs[result == "Positive", .SD[which.max(t)], by = id
][, .(me = signif(median(t), 2), 
      lo = signif(median(t) - IQR(t)/2, 2),
      hi = signif(median(t) + IQR(t)/2, 2)), by = VOC]

#--- extra figures with same data
# no of points available per patient
dt_ct_no_inf_before_start[, .N, by = c("data_id", "VOC", "result")][
  , VOC := fct_relevel(VOC, "Delta")] |> 
  ggplot() +
  geom_bar(aes(x = N, fill = VOC)) +
  theme_minimal() + 
  labs(x = "Number of PCR tests per id", y = "Count") + 
  facet_wrap(~result, scales = "free")

# delay between first viral load data and symptom onset
dt_ct_no_inf_before_start[, .(onset_time), by = c("data_id", "VOC")][
  , VOC := fct_relevel(VOC, "Delta")] |>
  unique() |> 
  ggplot() +
  geom_bar(aes(x = onset_time, fill = VOC)) +
  theme_minimal() + 
  labs(x = "Number of PCR tests per id", y = "Count")

# observed peak viral load
dt_ct_no_inf_before_start[, .SD[which.min(ct_value)], by = id
                          ][, VOC := fct_relevel(VOC, "Delta")] |>
  ggplot(aes(x = VOC, y = ct_value, colour = VOC, fill = VOC)) +
  geom_violin(alpha = 0.1) +
  geom_jitter(height = 0, width = 0.1) +
  theme_minimal() + 
  labs(x = "Number of PCR tests per id", y = "Count")

# time to clearance
dt_ct_no_inf_before_start[result == "Positive", .SD[which.max(t)], by = id
][, VOC := fct_relevel(VOC, "Delta")] |>
  ggplot(aes(x = VOC, y = t, colour = VOC, fill = VOC)) +
  geom_violin(alpha = 0.1) +
  geom_jitter(height = 0, width = 0.1) +
  theme_minimal() + 
  labs(x = "Number of PCR tests per id", y = "Count")


