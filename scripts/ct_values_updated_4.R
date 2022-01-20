library(data.table)
library(lubridate)
library(patchwork)
library(rstan)
library(cowplot)
library(stringr)

source("R/stan_data.R")
source("R/custom_plot_theme.R")
source("scripts/dry_vs_wet.R")

dt.ct <- fread("data/ct_values_updated_3.csv", drop = c("V7", "V8")) %>% 
  setnames(., c("Study number", "sample date", "swab type",
                "Ct", "VOC"), 
           c("id", "date", "swab_type", "ct_value",
             "voc")) %>% 
  .[, date := dmy(date)] %>% 
  .[voc == "Alpha", voc := "alpha"] %>% 
  .[voc == "Delta", voc := "delta"] %>% 
  .[voc == "Omicron", voc := "omicron"] %>% 
  .[voc == "OMicron", voc := "omicron"] %>% 
  .[voc == "Omicron/SGTF", voc := "omicron"] %>% 
  .[voc == "S gene present", voc := "omicron"] %>% 
  .[voc == "omicron" | voc == "delta" |
      voc == "Positive" | voc == "positive" |
      voc == "alpha" | voc == "Positive SARS-CoV-2", result := 1] %>% 
  .[voc == "negative" | voc == "SARS-CoV-2 Not Detected", result := 0] %>% 
  .[voc == "SARS-CoV-2 Not Detected" | 
      voc == "SARS-CoV-2 Inconclusive" | voc == "?SGTF" |
      voc == "retest/inconsistent" | voc == "Positive" | voc == "positive" |
      voc == "Positive SARS-CoV-2" | voc == "Inconclusive" |
      voc == "inconclusive", voc := NA] %>% 
  .[is.na(ct_value), ct_value := NA] %>% 
  .[result == 0, ct_value := 45] 

dt.symptoms <- fread("data/symptom_data_new.csv") %>% 
  setnames(., "start_date", "symptom_onset") %>% 
  .[, symptom_onset := dmy(symptom_onset)]

dt.all <- merge.data.table(dt.ct, dt.symptoms, all.x = TRUE,
                 by = c("id", "infection_id")) %>% 
  .[, t_first_test := as.numeric(date - min(date), units = "days"), 
    by = c("id", "infection_id")] %>% 
  .[, t_symp_onset := as.numeric(date - symptom_onset, units = "days"), 
    by = c("id", "infection_id")] %>% 
  .[, onset_rel :=  as.numeric(symptom_onset - min(date), units = "days"), 
    by = c("id", "infection_id")]

dt.all[swab_type == "VTM" & voc != "negative", ct_value_adjusted := ct_value*adjustment.dt[, mi]]
dt.all[swab_type != "VTM", ct_value_adjusted := ct_value]
dt.all[ct_value_adjusted == 45, ct_value_adjusted := NA]

dt.plot <- dt.all[,  `:=`  (result = factor(result,
                                           levels = c(1, 0),
                                           labels = c("Positive", "Negative")),
                           voc = factor(voc,
                                        levels = c("alpha", "delta", "omicron", "negative"),
                                        labels = c("Alpha", "Delta", "Omicron", "Negative")))]

#--- plots of dynamics since symptom onset - unadjusted
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
  .[t_symp_onset < 35 & t_symp_onset >= 0] %>%  
  ggplot(aes(x = t_symp_onset, y = ct_value)) + 
  geom_point(aes(colour = factor(voc), group = interaction(voc, id)), alpha = 0.4) + 
  geom_line(aes(colour = factor(voc),  group = interaction(voc, id)), alpha = 0) +
  geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
  stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since symptom onset - no adjustment", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1")

p.both <- p1 + p2 + plot_layout(guides = "collect")

ggsave("outputs/individual_ct_trajectories.png",
       p.both,
       width = 8,
       height = 4,
       bg = "white")

#--- plots of dynamics since symptom onset - adjusted
dt.plot.non.adj <- dt.plot[, c("t_symp_onset", "id", "voc", 
                               "ct_value", "t_first_test", "onset_rel", 
                               "infection_id")] %>% 
  .[, type := "not adjusted"]

dt.plot.adj <- dt.plot[swab_type == "VTM", ct_value := ct_value*adjustment.dt[, mi]] %>% 
  .[, c("t_symp_onset", "id", "voc", "ct_value", "t_first_test", 
        "onset_rel", "infection_id")] %>% 
  .[, type := "adjusted"]

dt.plot.both <- rbind(dt.plot.non.adj, dt.plot.adj) %>% 
  .[, type := factor(type, 
                     levels = c("adjusted", "not adjusted"), 
                     labels = c("Adjusted", "Not adjusted"))]

p_adj <- dt.plot.both[(voc == "Delta" | voc == "Omicron") &
                      type == "Adjusted" & is.na(ct_value) == FALSE] %>% 
  .[, obs_total := .N, by = "id"] %>% 
  #.[obs_total > 1] %>% 
  .[t_symp_onset < 30 & t_symp_onset >= 0] %>%  
  ggplot(aes(x = t_symp_onset, y = ct_value)) + 
  geom_point(aes(group = interaction(type, id)), alpha = 0.2) + 
  geom_line(aes(group = interaction(type, id)), alpha = 0) +
  geom_smooth(se = FALSE, ) + 
  stat_smooth(alpha = 0.1) +
  scale_y_reverse(lim = c(45, NA)) +
  coord_cartesian(xlim = c(0, 10)) +
  custom_plot_theme(flip = TRUE) +
  labs(title = "Since symptom onset", x = "Time (days)", y = "Ct value") +
  scale_colour_brewer(palette = "Set1") + 
  facet_wrap(~voc)

ggsave("outputs/individual_ct_trajectories_adjusted_new.png", 
       p_adj,
       height = 4,
       width = 8,
       bg = "white")

# p2_adj <- dt.plot.adj[(voc == "Alpha" | voc == "Delta" | voc == "Omicron") & is.na(ct_value) == FALSE] %>% 
#   .[, obs_total := .N, by = "id"] %>% 
#   #.[obs_total > 1] %>% 
#   .[t_symp_onset < 35 & t_symp_onset >= 0] %>%  
#   ggplot(aes(x = t_symp_onset, y = ct_value)) + 
#   geom_point(aes(colour = factor(voc), group = interaction(voc, id)), alpha = 0.4) + 
#   geom_line(aes(colour = factor(voc),  group = interaction(voc, id)), alpha = 0) +
#   geom_smooth(se = FALSE, aes(colour = factor(voc))) + 
#   stat_smooth(aes(colour = factor(voc)), alpha = 0.1) +
#   scale_y_reverse(lim = c(45, NA)) +
#   custom_plot_theme(flip = TRUE) +
#   labs(title = "Since symptom onset - adjusted ", x = "Time (days)", y = "Ct value") +
#   scale_colour_brewer(palette = "Set1")

p.both_adj <- p1_adj + p2_adj + plot_layout(guides = "collect")

ggsave("outputs/individual_ct_trajectories_adjusted.png",
       p.both,
       width = 8,
       height = 4,
       bg = "white")

# p3 <- dt.plot[(voc == "Delta" | voc == "Negative") & is.na(ct_value) == FALSE] %>% 
#   .[t_first_test >= 0] %>%  
#   .[t_symp_onset < 35 & t_symp_onset >= 0] %>%
#   .[, obs_total := .N, by = "id"] %>% 
#   .[obs_total > 1] %>% 
#   .[, .SD[!all(result == "Negative")], by = "id"] %>% 
#   ggplot(aes(x = t_first_test, y = ct_value)) + 
#   geom_point(aes(colour = factor(result), group = interaction(id, result)), alpha = 0.4) + 
#   geom_line(aes(colour = factor(result),  group = interaction(id, result)), alpha = 0.4) +
#   geom_vline(aes(xintercept = onset_rel), linetype = "dashed") +
#   #geom_smooth() +
#   facet_wrap(~id) +
#   custom_plot_theme(flip = TRUE) +
#   #theme(legend.position = "none") +
#   scale_y_reverse(lim = c(45, 10)) +
#   labs(title = "Delta", x = "Days since first swab", y = "Ct value") +
#   scale_color_manual(values = c("blue", "red"))
# 
# p4 <- dt.plot[(voc == "Omicron" | voc == "Negative") & is.na(ct_value) == FALSE] %>% 
#   .[t_first_test >= 0 & t_first_test <= 25] %>%  
#   .[t_symp_onset < 35 & t_symp_onset >= 0] %>%
#   .[, obs_total := .N, by = "id"] %>% 
#   .[obs_total > 1] %>% 
#   .[, .SD[!all(result == "Negative")], by = "id"] %>%  
#   ggplot(aes(x = t_first_test, y = ct_value)) + 
#   geom_point(aes(colour = factor(result), group = interaction(id, result)), alpha = 0.4) + 
#   geom_line(aes(colour = factor(result),  group = interaction(id, result)), alpha = 0.4) +
#   geom_vline(aes(xintercept = onset_rel), linetype = "dashed") + 
#   #geom_smooth() +
#   facet_wrap(~id) +
#   custom_plot_theme(flip = TRUE) +
#   #theme(legend.position = "none") +
#   scale_y_reverse(lim = c(45, 10)) +
#   labs(title = "Omicron", x = "Days since first swab", y = "Ct value") +
#   scale_color_manual(values = c("blue", "red"))
# 
# figure_2 <- p3 + p4 + plot_layout(guides = "collect")
# 
# ggsave("outputs/ct_values_individual_by_variant.png",
#        figure_2,
#        width = 15,
#        height = 7,
#        bg = "white")

#--- running inference

mn <- dt.all[, min(ct_value_adjusted, na.rm = TRUE)]
mx <- dt.all[, max(ct_value_adjusted, na.rm = TRUE)]

# scale Ct values
dt.all[, ct_scaled := (ct_value_adjusted - mn)/(mx - mn)]
#--- individual level data - NOT WORKING
# choosing with voc to fit to 
# dt.ct.delta <- dt.all[(voc == "delta" | voc == "negative") & is.na(ct_value) == FALSE] %>% 
#   .[t_first_test < 20 & t_first_test >= 0] %>% 
#   .[, obs_total := .N, by = c("id", "infection_id")] %>%
#   .[, .SD[obs_total > 2], by = id] %>% 
#   .[, day_rel := t_first_test + 5, units = "days", by = c("id", "infection_id")] %>%
#   .[, id := .GRP, by = id] %>% 
#   .[order(id, day_rel)] %>% 
#   .[order(id, day_rel, result),
#     c("id", "date", "ct_value", "ct_value_adjusted", "result", "day_rel")] 
# 
# dt.ct.omicron <- dt.all[(voc == "omicron" | voc == "negative") & is.na(ct_value) == FALSE] %>% 
#   .[t_first_test < 20 & t_first_test >= 0] %>% 
#   .[, obs_total := .N, by = c("id", "infection_id")] %>%
#   .[, .SD[obs_total > 2], by = id] %>% 
#   .[, day_rel := as.numeric(date - min(date) + 10, units = "days"), by = id] %>%
#   .[, id := .GRP, by = id] %>% 
#   .[order(id, day_rel)]

#--- pooled data
dt.ct.delta.pooled <- dt.all[(voc == "delta" | voc == "negative") &
                              t_first_test < 20 & t_first_test >= 0 &
                              is.na(ct_value) == FALSE] %>% 
  .[, id := .GRP, by = id] %>% 
  .[, c("id", "date", "infection_id", "ct_value",
        "ct_value_adjusted", "ct_scaled", "result", "t_first_test")] 

dt.ct.omicron.pooled <- dt.all[(voc == "omicron" | voc == "negative") &
                               t_first_test < 20 & t_first_test >= 0 &
                               is.na(ct_value) == FALSE] %>% 
  #.[, day_rel := as.numeric(date - min(date), units = "days"), by = id] %>%
  .[, id := .GRP, by = id] %>% 
  .[, c("id", "date", "infection_id", "ct_value",
        "ct_value_adjusted", "ct_scaled", "result", "t_first_test")] 

options(mc.cores = parallel::detectCores()) 
#--- choose here whether you want the individual-level fits 
#--- or a pooled fit
mod <- rstan::stan_model("stan/ct_trajectory_model_pooled.stan")

#--- fitting for delta
stan_data_delta <- stan_data_fun(dt.ct.delta.pooled)

fit_delta <- sampling(mod,
                      chains = 4,
                      data = stan_data_delta,
                      iter = 2000)

stan_trace(fit_delta, pars = c("t_p_mean", "t_s_mean", "t_lod_mean",
                               "c_p_mean", "c_s_mean", "sigma_obs",
                               "t_lod_abs"))

res_delta <- extract(fit_delta)

dt.c.peak.delta <- res_delta$c_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn]  %>% 
  .[, voc := "Delta"]

delta_c_peak <- dt.c.peak.delta %>% 
  ggplot() +
  geom_density(aes(smp), fill = "lightblue") + 
  #xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme()

dt.t.peak.delta <- res_delta$t_p_mean %>% 
  reshape2::melt() %>% 
  data.table()
  # .[, smp := (mx - mn) * value + mn]  %>% 
  # .[, voc := "Delta"]

delta_t_peak <- dt.t.peak.delta %>% 
  ggplot() +
  geom_density(aes(value), fill = "lightblue") + 
  #xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme()

extract(fit_delta, pars = "T_e") %>% 
  reshape2::melt(.) %>% 
  ggplot() +
  geom_density(aes(value)) + 
  xlim(0, 10) +
  labs(x = "Ct value", y = "Probability") + 
  facet_wrap(~Var2, scales = "free") + 
  custom_plot_theme()

# fit_dt_delta <- as.data.frame(fit_delta, pars = "ct") %>%
#   data.table() %>%
#   melt(., measure.vars = patterns("ct")) %>% 
#   .[, c("id", "time") := tstrsplit(variable, ",")] %>% 
#   .[, id := str_remove(id, "ct\\[")] %>% 
#   .[, time := str_remove(time, "]")] %>%
#   .[, iteration := 1:.N, by = "variable"] %>% 
#   .[, variable := NULL] %>% 
#   .[, time := as.numeric(time)] %>% 
#   .[, value := (mx - mn) * value + mn]

fit_dt_delta_pooled <- as.data.frame(fit_delta, pars = "ct") %>%
  data.table() %>%
  melt(., measure.vars = patterns("ct")) %>% 
  .[, time := sort(rep(1:20, 4000))] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)] %>%
  .[, value := (mx - mn) * value + mn]

fit_dt_summary <- fit_dt_delta[, .(me = quantile(value, 0.5), 
                                   lo = quantile(value, 0.025),
                                   hi = quantile(value, 0.975)),
                               by = c("id", "time")]

fit_dt_pooled_summary <- fit_dt_delta_pooled[, .(me = quantile(value, 0.5), 
                                                  lo = quantile(value, 0.025),
                                                  hi = quantile(value, 0.975)),
                                              by = time]

fit_dt_pooled_summary %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = me)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = "dodgerblue"), alpha = 0.1) +
  #facet_wrap(~id) + 
  scale_y_reverse(lim = c(30, NA)) +
  scale_y_reverse() +
  #xlim(0, 10) + 
  custom_plot_theme(flip = TRUE)

fit_delta@model_pars
stan_dens(fit_delta, pars = "c_p_mean")

#--- fitting for omicron
stan_data_omicron <- stan_data_fun(dt.ct.omicron.pooled)

fit_omicron <- sampling(mod,
                        chains = 4,
                        data = stan_data_omicron,
                        iter = 2000)

res_omicron <- extract(fit_omicron)

dt.peak.omicron <- res_omicron$c_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn] %>% 
  .[, voc := "Omicron"]

omicron_c_peak <- dt.peak.omicron %>% 
  ggplot() +
  geom_density(aes(smp), fill = "lightblue") + 
  xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme()

rbind(dt.peak.delta, dt.peak.omicron) %>% 
  ggplot() +
  geom_density(aes(smp, fill = voc)) + 
  xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme()

extract(fit_omicron, pars = "T_e") %>% 
  reshape2::melt(.) %>% 
  ggplot() +
  geom_density(aes(value,)) + 
  xlim(0, 6) +
  labs(x = "Ct value", y = "Probability") + 
  facet_wrap(~Var2, scales = "free") + 
  custom_plot_theme()

hello <- as.data.frame(fit_omicron, pars = "ct") %>%
  data.table() %>%
  melt(., measure.vars = patterns("ct")) %>% 
  .[, time := sort(rep(1:20, 1000))] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)] %>%
  .[, value := (mx - mn) * value + mn]

fit_dt_omicron_pooled_summary <- fit_dt_omicron_pooled[, .(me = quantile(value, 0.5), 
                                                           lo = quantile(value, 0.025),
                                                           hi = quantile(value, 0.975)),
                                                       by = time]

p2 <- fit_dt_omicron_pooled_summary %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = me)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = "dodgerblue"), alpha = 0.1) +
  #facet_wrap(~id) + 
  scale_y_reverse(lim = c(30, NA)) +
  scale_y_reverse() +
  #xlim(0, 10) + 
  custom_plot_theme(flip = TRUE)

#---

ggsave("delta_c_peak.png",
       delta_c_peak,
       width = 6,
       height = 6,
       bg = "white")

delta_t_peak <- dt.tmp %>% 
  ggplot() +
  geom_density(aes(smp), fill = "lightblue") + 
  xlim(0, 45) +
  labs(x = "Time (days)", y = "Probability",
       title = "Posterior distribution for timing of peak Ct value - pooled for all individuals") +
  custom_plot_theme()


stan_dens(fit_delta, pars = "c_p_mean") + 
  xlim(0, 5)

res <- extract(fit_delta)

# allsamp <- data.table::melt(as.data.table(res$T_e)[, iterations := 1:.N],
#                             id.vars = c("iterations"),
#                             value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
#                             ][, num_id := as.numeric(num_id)]
# 
# 
# tmp_plot <- allsamp %>%
#   ggplot() +
#   geom_density(aes(x = smp)) +
#   facet_wrap(~num_id)

# tmp <- data.table::melt(as.data.table(res_$c_p)[, iterations := 1:.N],
#                         id.vars = c("iterations"),
#                         value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
#                         ][, num_id := as.numeric(num_id)] %>%
#   .[, smp := (mx - mn) * smp + mn]

fit_dt <- as.data.frame(fit_omicron, pars = "ct") %>%
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
  #geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = "dodgerblue"), alpha = 0.1) +
  facet_wrap(~id) + 
  scale_y_reverse(lim = c(30, NA)) +
  xlim(0, 10) + 
  custom_plot_theme(flip = TRUE)

tmp

