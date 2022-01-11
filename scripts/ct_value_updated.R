library(data.table)
library(lubridate)
library(patchwork)
library(rstan)

source("ct_dynamics/R/stan_data.R")
source("ct_dynamics/R/stan_data.R")

dt.ct <- fread("ct_dynamics/data/ct_values_updated.csv") %>% 
  setnames(., c("Study number", "sample date", "barcode", 
                "Ct", "swab type", "sample order", "Day post swab 1",
                "VOC", "Symptom onset", "Days post symptom onset"), 
           c("id", "date", "barcode", "ct_value", "swab_type", "sample_id",
             "days_since_first_swab", "voc", "symptom_onset", 
             "symptom_onset_rel")) %>% 
  .[, date := dmy(date)] %>% 
  .[, symptom_onset := dmy(symptom_onset)] %>% 
  .[voc == "Delta", voc := "delta"] %>% 
  .[voc == "Omicron", voc := "omicron"] %>% 
  .[voc == "OMicron", voc := "omicron"] %>% 
  .[voc == "omicron" | voc == "delta" |
    voc == "Positive" | voc == "positive" |
    voc == "alpha" | voc == "Positive SARS-CoV-2", result := 1] %>% 
  .[voc == "negative" | voc == "SARS-CoV-2 Not Detected", result := 0] %>% 
  .[voc == "SARS-CoV-2 Not Detected" | 
    voc == "SARS-CoV-2 Inconclusive" | voc == "?SGTF" |
    voc == "retest/inconsistent" | voc == "Positive" | voc == "positive" |
    voc == "Positive SARS-CoV-2", voc := NA] %>% 
  .[is.na(ct_value), ct_value := NA] %>% 
  .[result == 0, ct_value := 45] %>% 
  .[, day_rel := as.numeric(date - min(date), units = "days"), 
    by = c("id", "infection_id")] %>% 
  .[, symp_rel := as.numeric(symptom_onset - min(date), units = "days"), 
    by = c("id", "infection_id")]

# dt.symptom.onset <- fread("ct_dynamics/data/symptom_onset.csv") %>% 
#   .[, symptom_onset := dmy(symptom_onset)] %>% 
#   na.omit(.) %>% 
#   setnames(., c("participant_id", "infection_number"), c("id", "infection_id"))

dt.plot <- dt.all[, result := factor(result,
                                     levels = c(1, 0), 
                                     labels = c("Positive", "Negative"))]

# dt.all[,`:=` (day_rel = day_rel + abs(symp_rel), symp_rel = symp_rel - symp_rel), 
#        by = c("id", "infection_id")]

p1 <- dt.plot[(voc == "delta" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[day_rel < 25] %>%  
  ggplot(aes(x = day_rel, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = interaction(id, result)), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = interaction(id, result)), alpha = 0.4) +
  geom_smooth() +
  #facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  #scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

p2 <- dt.plot[(voc == "omicron" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
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
  #scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Omicron", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

figure_1 <- p1 + p2 + plot_layout(guides = "collect")

ggsave("ct_dynamics/outputs/ct_dynamics_pooled_by_variant.png",
       figure_1,
       width = 8,
       height = 4,
       bg = "white")

p3 <- dt.plot[(voc == "delta" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[days_since_first_swab < 25] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = days_since_first_swab, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  geom_vline(aes(xintercept = symp_rel)) +
  #geom_smooth() +
  facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

p4 <- dt.plot[(voc == "omicron" | is.na(voc) == TRUE) & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  .[day_rel < 12] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = day_rel, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  geom_vline(aes(xintercept = symp_rel)) + 
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

# preparing stan data
mn <- dt.ct[, min(ct_value, na.rm = TRUE)]
mx <- dt.ct[, max(ct_value, na.rm = TRUE)]

# scale Ct values
dt.ct[, ct_scaled := (ct_value - mn)/(mx - mn)]

# choosing with voc to fit to 
dt.ct.delta <- dt.ct[voc == "delta" | voc == "negative"] %>%
  .[, id := .GRP, by = id]

start_date_delta <- dt.ct.delta[, min(symptom_onset, na.rm = TRUE) - 1]

dt.ct.delta <- dt.ct.delta[is.na(result) == FALSE] %>% 
  .[, day_rel := date - start_date_delta, by = "id"] %>% 
  .[is.na(symptom_onset), symptom_onset := min(date) - 7,
    by = "id"] %>% 
  .[, symp_rel := as.numeric(symptom_onset - start_date_delta, units = "days"),
    by = "id"] %>% 
  .[, te_upper_bound := ifelse(
    any(day_rel[result == 1] < symp_rel[result == 1]),
    min(day_rel[result == 1 & day_rel < symp_rel]),
    symp_rel), by = id] %>% 
  .[, te_upper_bound := unique(te_upper_bound), id]

# Construct stan data
stan_data_delta <- stan_data_fun(dt.ct.delta)
# stan_data_omicron <- stan_data_fun(dt.ct.omicron)
# stan_data_both <- stan_data_fun(dt.together_both)

options(mc.cores = parallel::detectCores()) 
# rstan_options(auto_write = TRUE)
mod <- rstan::stan_model("ct_dynamics/stan/ct_trajectory_model.stan")

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

tmp <- stan_dens(fit_delta, pars = "T_e") 

ggsave("ct_dynamics/outputs/all_Te.png", 
       tmp,
       width = 10,
       height = 10,
       bg = "white")
