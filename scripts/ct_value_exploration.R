library(data.table)
library(lubridate)
library(patchwork)
library(rstan)

source("ct_dynamics/R/stan_data.R")

dt.ct <- fread("ct_dynamics/data/ct_value_by_variant.csv") %>% 
  setnames(., c("Study number", "sample date", "barcode", 
                "Ct", "sample order", "Day post swab 1", "VOC"), 
           c("id", "date", "barcode", "ct_value", "sample_id",
             "days_since_first_swab", "voc")) %>% 
  .[, date := dmy(date)] %>% 
  .[voc == "Delta", voc := "delta"] %>% 
  .[voc == "Omicron", voc := "omicron"] %>% 
  .[voc == "OMicron", voc := "omicron"] %>% 
  .[voc == "omicron" | voc == "delta", result := 1] %>% 
  .[voc == "negative", result := 0] %>% 
  .[voc == "negative", voc := NA] %>% 
  .[is.na(ct_value), ct_value := NA]

dt.symptom.onset <- fread("ct_dynamics/data/symptom_onset.csv") %>% 
  .[, symptom_onset := dmy(symptom_onset)] %>% 
  na.omit(.) %>% 
  setnames(., c("participant_id", "infection_number"), c("id", "infection_id"))

dt.all <- merge.data.table(dt.ct, dt.symptom.onset, all.x = TRUE, by = c("id", "infection_id")) %>% 
  .[result == 0, ct_value := 45] %>% 
  .[, result := factor(result,
                       levels = c(1, 0), 
                       labels = c("Positive", "Negative"))] %>% 
  .[, day_rel := as.numeric(date - min(date), units = "days"), by = c("id", "infection_id")] %>% 
  .[, symp_rel := as.numeric(symptom_onset - min(date), units = "days"), by = c("id", "infection_id")]


dt.all[,`:=` (day_rel = day_rel + abs(symp_rel), symp_rel = symp_rel - symp_rel), 
       by = c("id", "infection_id")]

p1 <- dt.all[(voc == "delta" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[days_since_first_swab < 25] %>%  
  ggplot(aes(x = days_since_first_swab, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = interaction(id, result)), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = interaction(id, result)), alpha = 0.4) +
  geom_smooth() +
  #facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

p2 <- dt.all[(voc == "omicron" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[days_since_first_swab < 12] %>%  
  ggplot(aes(x = days_since_first_swab, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = interaction(id, result)), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = interaction(id, result)), alpha = 0.4) +
  geom_smooth() +
  #facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Omicron", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

figure_1 <- p1 + p2 + plot_layout(guides = "collect")

ggsave("ct_dynamics/outputs/ct_dynamics_pooled_by_variant.png",
       figure_1,
       width = 8,
       height = 4,
       bg = "white")

p3 <- dt.all[(voc == "delta" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[days_since_first_swab < 25] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = days_since_first_swab, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  # geom_vline(aes(xintercept = symp_rel))
  #geom_smooth() +
  facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

p4 <- dt.all[(voc == "omicron" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[days_since_first_swab < 12] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = days_since_first_swab, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  #geom_smooth() +
  facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Omicron", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

figure_2 <- p3 + p4 + plot_layout(guides = "collect")

ggsave("ct_dynamics/outputs/ct_dynamics_individual_by_variant.png",
       figure_2,
       width = 12,
       height = 7,
       bg = "white")

dt.together <- merge.data.table(dt.ct, dt.symptom.onset, by = c("participant_id"))
dt.together[participant_id == 883 & date == ymd("2021-06-14"), symptom_onset := ymd("2021-06-10")]
dt.together[sample_id == "RLNA32883B", c("voc", "day_rel") := list("omicron", 6)]
dt.together[, symptom_onset_rel := as.numeric(symptom_onset - min(date), units = "days"),
            by = c("participant_id", "voc")]

# imputing missing onset dates where none are given
dt.together[is.na(symptom_onset_rel) == TRUE, symptom_onset_rel := -7]

# dt.together[is.na(symptom_onset_rel) == TRUE, 
#             symptom_onset_rel := mean(unique(symptom_onset_rel), na.rm = TRUE)]

# preparing stan data
mn <- dt.together[, min(ct_value, na.rm = TRUE)]
mx <- dt.together[, max(ct_value, na.rm = TRUE)]

# scale Ct values
dt.together[, ct_scaled := (ct_value - mn)/(mx - mn)]

dt.together[voc == "delta" | voc == "omicron", result := 1]
dt.together[voc == "negative", result := 0]
dt.together <- dt.together[voc != "? SGTF Sept"]
#dt.together <- dt.together[day_rel != 189]

# dt.together <- dt.together[obs > 1]

# choosing with voc to fit to 
dt.together_delta <- dt.together[voc == "delta" | voc == "negative"] %>%
  .[, id := .GRP, by = participant_id] 

dt.together_omicron <- dt.together[voc == "omicron" | voc == "negative"] %>%
  .[, id := .GRP, by = participant_id]

dt.together_both <- dt.together[, id := .GRP, by = participant_id]



# Construct stan data
stan_data_delta <- stan_list(dt.together_delta)
stan_data_omicron <- stan_list(dt.together_omicron)
stan_data_both <- stan_list(dt.together_both)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mod <- stan_model("stan/ct_trajectory_model.stan")

fit_delta <- sampling(mod,
                      data = stan_data_delta,
                      iter = 2000)

fit_omicron <- sampling(mod,
                        data = stan_data_omicron,
                        iter = 2000)

fit_both <- sampling(mod,
                     data = stan_data_both,
                     iter = 4000)

res_both <- extract(fit_both)

tmp <- rstan::traceplot(fit_omicron, pars = "T_e")

ggsave("trace_plots_delta.png", 
       tmp,
       width = 10,
       height = 10,
       bg = "white")

stan_dens(fit_delta, pars = "c_s")
stan_dens(fit_omicron, pars = "t_lod")
stan_dens(fit_both, pars = "t_lod")

tmp <- stan_dens(fit_both, pars = "T_e")
ggsave("T_e_plots.png", 
       tmp,
       width = 10,
       height = 10,
       bg = "white")
