library(data.table)
library(lubridate)
library(patchwork)
library(rstan)
library(cowplot)
library(stringr)

source("R/stan_data.R")
# source("R/stan_data.R")
source("R/custom_plot_theme.R")

# dt.ct <- fread("data/ct_values_updated_2.csv") %>% 
#   setnames(., c("Study number", "sample date", 
#                 "Ct", "VOC", "symptom onset check EW"), 
#            c("id", "date", "ct_value", "voc", "symptom_onset")) %>% 
#   .[, date := dmy(date)] %>% 
#   .[, symptom_onset := dmy(symptom_onset)] %>% 
#   .[voc == "Delta", voc := "delta"] %>% 
#   .[voc == "Omicron", voc := "omicron"] %>% 
#   .[voc == "OMicron", voc := "omicron"] %>% 
#   .[voc == "omicron" | voc == "delta" |
#       voc == "Positive" | voc == "positive" |
#       voc == "alpha" | voc == "Positive SARS-CoV-2", result := 1] %>% 
#   .[voc == "negative" | voc == "SARS-CoV-2 Not Detected", result := 0] %>% 
#   .[voc == "SARS-CoV-2 Not Detected" | 
#     voc == "SARS-CoV-2 Inconclusive" | voc == "?SGTF" |
#     voc == "retest/inconsistent" | voc == "Positive" | voc == "positive" |
#     voc == "Positive SARS-CoV-2", voc := NA] %>% 
#   .[is.na(ct_value), ct_value := NA] %>% 
#   .[result == 0, ct_value := 45] %>% 
#   .[, days := as.numeric(date - min(date), units = "days"), 
#     by = c("id", "infection_id")] %>% 
#   .[, symp_rel := as.numeric(symptom_onset - min(date), units = "days"), 
#     by = c("id", "infection_id")]

dt.ct <- fread("data/ct_values_updated_2.csv") %>% 
  setnames(., c("Study number", "sample date", 
                "Ct", "VOC", "symptom onset check EW"), 
           c("id", "date", "ct_value", "voc", "symptom_onset")) %>% 
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
  .[, t_first_test := as.numeric(date - min(date), units = "days"), 
    by = c("id", "infection_id")] %>% 
  .[, t_symp_onset := as.numeric(date - symptom_onset, units = "days"), 
    by = c("id", "infection_id")]

# dt.ct.tmp.1 <- dt.ct[, .SD[symptom_onset < min(date)], by = c("id", "infection_id")] %>%
#   .[, `:=` (day_rel  = min(t_first_test) - t_symp_onset + 2,
#             symp_rel = t_symp_onset - t_symp_onset + 2), by = c("id", "infection_id", "date")]
# 
# dt.ct.tmp.2 <- dt.ct[, .SD[symptom_onset >= min(date) | is.na(symptom_onset) == TRUE],
#                      by = c("id", "infection_id")] %>%
#   .[, `:=` (day_rel  = min(t_first_test) + 2,
#             symp_rel = t_symp_onset + 2), by = c("id", "infection_id", "date")]
# 
# dt.ct <- rbind(dt.ct.tmp.1, dt.ct.tmp.2)[order(id, date, voc)]

dt.plot <- dt.ct[,  `:=`  (result = factor(result,
                                          levels = c(1, 0), 
                                          labels = c("Positive", "Negative")),
                           voc = factor(voc,
                                        levels = c("alpha", "delta", "omicron", "negative"),
                                        labels = c("Alpha", "Delta", "Omicron", "Negative")))]

p1 <- dt.plot[(voc == "Alpha" | voc == "Delta" | voc == "Omicron") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  #.[obs_total > 1] %>% 
  .[t_first_test < 30] %>%  
  ggplot(aes(x = t_first_test, y = ct_value)) + 
  geom_point(aes(colour = factor(voc), group = interaction(id, voc)), alpha = 0.4) + 
  geom_line(aes(colour = factor(voc),  group = interaction(id, voc)), alpha = 0) +
  geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
  stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since first test", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1")

p2 <- dt.plot[(voc == "Alpha" | voc == "Delta" | voc == "Omicron") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  #.[obs_total > 1] %>% 
  .[t_symp_onset < 35] %>%  
  ggplot(aes(x = t_symp_onset, y = ct_value)) + 
  geom_point(aes(colour = factor(voc), group = interaction(voc, id)), alpha = 0.4) + 
  geom_line(aes(colour = factor(voc),  group = interaction(voc, id)), alpha = 0) +
  geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
  stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since symptom onset", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1")

figure_1 <- p1 + p2 + plot_layout(guides = "collect")
figure_1

p1 <- dt.plot[(voc == "Alpha" | voc == "Delta" | voc == "Omicron") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  #.[obs_total > 1] %>% 
  .[t_first_test < 30] %>%  
  ggplot(aes(x = t_first_test, y = ct_value)) + 
  geom_point(aes(colour = factor(voc), group = interaction(id, voc)), alpha = 0.4) + 
  geom_line(aes(colour = factor(voc),  group = interaction(id, voc)), alpha = 0) +
  geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
  stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since first test", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1")

p2 <- dt.plot[(voc == "Alpha" | voc == "Delta" | voc == "Omicron") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  #.[obs_total > 1] %>% 
  .[t_symp_onset < 35] %>%  
  ggplot(aes(x = t_symp_onset, y = ct_value)) + 
  geom_point(aes(colour = factor(voc), group = interaction(voc, id)), alpha = 0.4) + 
  geom_line(aes(colour = factor(voc),  group = interaction(voc, id)), alpha = 0) +
  geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
  stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since symptom onset", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1")

figure_1 <- p1 + p2 + plot_layout(guides = "collect")
figure_1

ggsave("outputs/ct_dynamics.png",
       figure_1,
       width = 8,
       height = 4,
       bg = "white")

p3 <- dt.plot[(voc == "delta" | voc == "negative") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  #.[day_rel < 25] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = day_rel_2, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  geom_vline(aes(xintercept = symp_rel_2), linetype = "dashed") +
  #geom_smooth() +
  facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

p4 <- dt.plot[(voc == "omicron" | voc == "negative") & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  .[obs_total > 1] %>% 
  # .[day_rel < 12] %>%  
  .[, .SD[!all(result == "Negative")], by = "id"] %>% 
  ggplot(aes(x = day_rel_2, y = ct_value)) + 
  geom_point(aes(colour = factor(result), group = id), alpha = 0.4) + 
  geom_line(aes(colour = factor(result),  group = id), alpha = 0.4) +
  geom_vline(aes(xintercept = symp_rel_2), linetype = "dashed") + 
  #geom_smooth() +
  facet_wrap(~id) +
  custom_plot_theme +
  #theme(legend.position = "none") +
  scale_y_reverse(lim = c(45, 10)) +
  labs(title = "Omicron", x = "Days since first swab", y = "Ct value") +
  scale_color_manual(values = c("blue", "red"))

figure_2 <- p3 + p4 + plot_layout(guides = "collect")

ggsave("outputs/ct_values_individual_by_variant.png",
       figure_2,
       width = 12,
       height = 7,
       bg = "white")

mn <- dt.ct[, min(ct_value, na.rm = TRUE)]
mx <- dt.ct[, max(ct_value, na.rm = TRUE)]

# scale Ct values
dt.ct[, ct_scaled := (ct_value - mn)/(mx - mn)]

# just for now, imputing most common symptom onset date,
# where none is present - need to think about how to deal with these 
# individuals properly down the line
dt.ct[is.na(symp_rel_2) == TRUE, symp_rel_2 := 2]

# choosing with voc to fit to 
dt.ct.delta <- dt.ct[voc == "delta" | voc == "negative"] %>%
  .[, id := .GRP, by = id] %>% 
  .[, `:=` (day_rel_2 = day_rel_2 + 10, 
            symp_rel_2 = symp_rel_2 + 10)]

stan_data_delta <- stan_data_fun(dt.ct.delta)

options(mc.cores = parallel::detectCores()) 
mod <- rstan::stan_model("stan/ct_trajectory_model.stan")

fit_delta <- sampling(mod,
                      chains = 1,
                      data = stan_data_delta,
                      iter = 2000)

res <- extract(fit_delta)

allsamp <- data.table::melt(as.data.table(res$T_e)[, iterations := 1:.N], 
                            id.vars = c("iterations"), 
                            value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
                            ][, num_id := as.numeric(num_id)]

tmp <- data.table::melt(as.data.table(res$c_p)[, iterations := 1:.N], 
                 id.vars = c("iterations"), 
                 value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
                 ][, num_id := as.numeric(num_id)] %>% 
  .[, smp := (mx - mn) * smp + mn]


tmpplot <- tmp %>% 
  ggplot() + 
  geom_density(aes(x = smp)) +
  facet_wrap(~num_id, scales = "free")

ggsave("outputs/all_cp_plots.png",
       tmpplot,
       width = 8,
       height = 10)


tmp <- stan_dens(fit_delta, pars = "c_p") + 
  xlim(0, 1)

res$c_p %>% melt(., measure.vars = patterns("V"))

fit_delta@model_pars

fit_dt <- as.data.frame(fit_delta, pars = "ct") %>%
  data.table() %>%
  melt(., measure.vars = patterns("ct")) %>% 
  .[, c("id", "time") := tstrsplit(variable, ",")] %>% 
  .[, id := str_remove(id, "ct\\[")] %>% 
  .[, time := str_remove(time, "]")] %>%
  .[, iteration := 1:.N, by = "variable"] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)] %>% 
  .[, value := (mx - mn) * value + mn]
 
fit_dt_summary <- fit_dt[, .(me = quantile(value, 0.5), 
                             lo = quantile(value, 0.025),
                             hi = quantile(value, 0.975)),
                         by = c("id", "time")]

tmp <- fit_dt_summary %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = lo)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = "dodgerblue"), alpha = 0.1) +
  facet_wrap(~id)

ggsave("outputs/all_cpeak_fits.png",
       tmp,
       width = 10,
       height = 10,
       bg = "white")
