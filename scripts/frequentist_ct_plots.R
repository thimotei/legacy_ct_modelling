library(ggsci)

dt_raw <- fread("data_raw/raw_data.csv")

dt_proc <- process_data(dt_raw)
dt_proc <- voc_status_attribution(dt_proc)
dt_proc <- t_last_exposure(dt_proc, imm_delay = 14)

obs <- subset_data(dt_proc,
                   no_pos_swabs = 2,
                   first_pos_min = -15)

obs[, unique(VOC)]

cols_to_remove <- c("centre", "sex", "age", "barcode", "Sequenced",
                    "dose_1", "dose_2", "dose_3", "dose_4")

cols_to_keep <- setdiff(colnames(obs), cols_to_remove)

dt_clean <- obs[, ..cols_to_keep] 

dt_plot <- obs[result == "Positive" & !VOC == "unknown"]




dt_no_points <- obs[, .N, by = c("id", "VOC", "symptoms", "age_group", "no_exposures")]

# dt_no_points[VOC == "Omicron (BA.1)", VOC := "BA.1"]
# dt_no_points[VOC == "Omicron (BA.2)", VOC := "BA.2"] 
# dt_no_points[symptoms == "symptomatic", symptoms := "+"]
# dt_no_points[symptoms == "asymptomatic", symptoms := "-"]

dt_no_points[, VOC := fct_relevel(VOC, "Delta")]
dt_no_points[, symptoms := fct_relevel(symptoms, "+")]
dt_no_points[, age_group := fct_relevel(age_group, c("20-34", "35-49", "50+"))]
dt_no_points[, no_exposures := fct_relevel(no_exposures, "3", "4", "5+")]

p1 <- dt_no_points |> 
  ggplot() +
  geom_bar(aes(x = N, fill = VOC)) +
  theme(legend.position = "bottom") +
  labs(fill = "")

p2 <- dt_no_points |> 
  ggplot() +
  geom_bar(aes(x = N, fill = symptoms)) +
  theme(legend.position = "bottom") +
  labs(fill = "")

p3 <- dt_no_points |> 
  ggplot() +
  geom_bar(aes(x = N, fill = age_group)) +
  theme(legend.position = "bottom") +
  labs(fill = "")

p4 <- dt_no_points |> 
  ggplot() +
  geom_bar(aes(x = N, fill = no_exposures)) +
  theme(legend.position = "bottom") +
  labs(fill = "")

p4
