library(data.table)
library(lubridate)
library(patchwork)
library(rstan)
library(cowplot)
library(stringr)

source("R/stan_data.R")
source("R/custom_plot_theme.R")
source("scripts/dry_vs_wet.R")

dt.raw <- fread("data/ct_values_clean.csv") 

dt.ct <- fread("data/ct_values_clean.csv") %>% 
  .[, swab_date := dmy(swab_date)] %>% 
  .[barcode %like% "49U", swab_date := swab_date - 1] %>% 
  .[symptom_onset_date == "unknown", symptom_onset_date := NA] %>% 
  .[, symptom_onset_date := dmy(symptom_onset_date)] %>% 
  .[, date_dose_1 := dmy(date_dose_1)] %>% 
  .[, date_dose_2 := dmy(date_dose_2)] %>% 
  .[, date_dose_3 := dmy(date_dose_3)] %>% 
  .[, days_since_symptom_onset := as.numeric(swab_date - symptom_onset_date, units = "days"), by = "ID"] %>% 
  .[, days_since_first_test := as.numeric(swab_date - min(swab_date), units = "days"), by = "ID"] %>% 
  .[ct != "unknown"] %>% 
  .[, ct := as.numeric(ct)] %>% 
  .[result == "Negative", ct := 40] %>% 
  .[swab_type == "VTM" & result != "Negative", ct_adjusted := ct*adjustment.dt[, mi]] %>%
  .[swab_type != "VTM", ct_adjusted := ct] %>% 
  .[result == "Positive" | result == "Inconclusive", result_num := 1] %>% 
  .[result == "Negative" | result == "Invalid" | result == "", result_num := 0]

no_pos_cts <- dt.ct[(result == "Positive" | result == "Inconclusive") & is.na(ct) == FALSE] %>% 
  .[, .N, by = c("ID")]

dt.ct <- merge.data.table(dt.ct, no_pos_cts, by = "ID") %>%  
  setnames(., 
           c("N", "number_vaccines (14 days pre ix)"), 
           c("no_pos_results", "no_vaccines"))

#--- removing duplicate VTM vs Dry swabs
barcodes.to.remove <-  dt.ct[swab_type == "Dry" | swab_type == "VTM"] %>% 
  .[, total_obs_by_date := .N, by = c("ID", "swab_date")] %>% 
  .[total_obs_by_date > 1, c("ID", "swab_date", "barcode", "swab_type")] %>% 
  .[order(ID, swab_date)] %>% 
  .[swab_type == "Dry", barcode]

dt.ct <- dt.ct[!barcode %in% barcodes.to.remove] 

#--- conditioning on only tests performed in the Crick, leaving out
#--- UCLH swabs for now
dt.ct <- dt.ct[centre == "crick"]
#--- running inference

mn <- dt.ct[, min(ct_adjusted, na.rm = TRUE)]
mx <- dt.ct[, max(ct_adjusted, na.rm = TRUE)]

# scale Ct values
dt.ct[, ct_scaled := (ct_adjusted - mn)/(mx - mn)]

#--- pooled data
dt.delta.pooled <- dt.ct[VOC == "Delta" & 
                        days_since_first_test < 20 & 
                        days_since_first_test >= 0] %>% 
  .[, ID := .GRP, by = ID]
  # .[, days_since_first_test := days_since_first_test + 15, by = ID]
 
# dt.delta.pooled[, c("ID", "ct", "days_since_first_test", "result")] %>%
#   .[(result == "Positive" | result == "Inconclusive") & is.na(ct) == FALSE,
#     no_pos_results := .N, by = ID] %>%
#   .[(result == "Positive" | result == "Inconclusive")] %>%
#   .[, .SD[no_pos_results > 2], by = ID] %>%
#   ggplot(aes(x = days_since_first_test, y = ct)) +
#   geom_point() +
#   geom_smooth()
#   facet_wrap(~ID)

dt.omicron.pooled <- dt.ct[VOC == "Omicron" & 
                        days_since_first_test < 20 & 
                        days_since_first_test >= 0] %>% 
  .[, ID := .GRP, by = ID]
  # .[, days_since_first_test := days_since_first_test, by = ID]

options(mc.cores = parallel::detectCores()) 
#--- choose here whether you want the individual-level fits 
#--- or a pooled fit
mod <- rstan::stan_model("stan/ct_trajectory_model_pooled.stan")

#--- fitting for delta
stan_data_delta <- stan_data_fun(dt.delta.pooled)

fit_delta <- sampling(mod,
                      chains = 4,
                      data = stan_data_delta,
                      iter = 2000)

pairs(fit_delta, pars = c("t_p_mean", "t_s_mean", "t_lod_mean",
                          "c_p_mean", "c_s_mean", "sigma_obs","t_lod_abs"))

# stan_trace(fit_delta, pars = c("t_p_mean", "t_s_mean", "t_lod_mean",
#                                "c_p_mean", "c_s_mean", "sigma_obs",
#                                "t_lod_abs"))

res_delta <- extract(fit_delta)

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

#--- fitting for omicron
stan_data_omicron <- stan_data_fun(dt.omicron.pooled)

fit_omicron <- sampling(mod,
                        chains = 4,
                        data = stan_data_omicron,
                        iter = 2000)

res_omicron <- extract(fit_omicron)

#--- combining fits
#--- delta fits
dt.t.e.delta <- res_delta$T_e %>% 
  reshape2::melt() %>% 
  data.table()  %>% 
  .[, voc := "Delta"]

dt.c.peak.delta <- res_delta$c_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn]  %>% 
  .[, voc := "Delta"]

dt.t.peak.delta <- res_delta$t_p_mean %>% 
  reshape2::melt() %>% 
  data.table()  %>% 
  .[, voc := "Delta"]

dt.t.e.omicron <- res_omicron$T_e %>% 
  reshape2::melt() %>% 
  data.table()  %>% 
  .[, voc := "Omicron"]

dt.c.peak.omicron <- res_omicron$c_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn] %>%
  .[, voc := "Omicron"]

dt.t.peak.omicron <- res_omicron$t_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, voc := "Omicron"]

T.e.posteriors <- rbind(dt.t.e.delta, dt.t.e.omicron) %>%
  ggplot() +
  geom_density(aes(value, fill = voc), alpha = 0.3) +
  # xlim(0, 5) +
  labs(x = "Time", y = "Probability",
       title = "Posterior distribution for exposure time") +
  custom_plot_theme() +
  scale_fill_brewer(palette = "Set2")

p.c.peak.posteriors <- rbind(dt.c.peak.delta, dt.c.peak.omicron) %>% 
  ggplot() +
  geom_density(aes(smp, fill = voc), alpha = 0.3) + 
  xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme() +
  scale_fill_brewer(palette = "Set2")

p.t.peak.posteriors <- rbind(dt.t.peak.delta, dt.t.peak.omicron) %>% 
  ggplot() +
  geom_density(aes(value, fill = voc), alpha = 0.3) + 
  xlim(0, 10) +
  labs(x = "Time (days after exposure)", y = "Probability",
       title = "Posterior distribution for timing of peak") +
  custom_plot_theme() +
  scale_fill_brewer(palette = "Set2")

p.all.posteriors <- p.c.peak.posteriors + p.t.peak.posteriors + plot_layout(guides = 'collect')

#--- plotting posterior predictive
fit_dt_delta_pooled <- as.data.frame(fit_delta, pars = "ct") %>%
  data.table() %>%
  melt(., measure.vars = patterns("ct")) %>% 
  .[, time := sort(rep(1:201, 4000))] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)] %>%
  .[, value := (mx - mn) * value + mn]

fit_dt_delta_pooled_summary <- fit_dt_delta_pooled[, .(me = quantile(value, 0.5), 
                                                       lo = quantile(value, 0.025),
                                                       hi = quantile(value, 0.975)),
                                                   by = time] %>% 
  .[, voc := "Delta"] %>% 
  .[, time := time/10 - 1]

fit_dt_omicron_pooled <- as.data.frame(fit_omicron, pars = "ct") %>%
  data.table() %>%
  melt(., measure.vars = patterns("ct")) %>% 
  .[, time := sort(rep(1:201, 4000))] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)] %>%
  .[, value := (mx - mn) * value + mn]

fit_dt_omicron_pooled_summary <- fit_dt_omicron_pooled[, .(me = quantile(value, 0.5), 
                                                           lo = quantile(value, 0.025),
                                                           hi = quantile(value, 0.975)),
                                                       by = time] %>% 
  .[, voc := "Omicron"] %>% 
  .[, time := time/10 - 1]

# 
# dt.ct.delta.pooled.plot <- dt.ct.delta.pooled[, ct_value_plot := ct_value_adjusted] %>% 
#   .[is.na(ct_value_adjusted), ct_value_adjusted := 45]
# 
# dt.ct.omicron.pooled.plot <- dt.ct.omicron.pooled[, ct_value_plot := ct_value_adjusted] %>% 
#   .[is.na(ct_value_adjusted), ct_value_adjusted := 45]

p.predictive <- rbind(fit_dt_delta_pooled_summary, fit_dt_omicron_pooled_summary) %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = me, colour = voc)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = voc), alpha = 0.2) +
  scale_y_reverse() +
  custom_plot_theme(flip = TRUE) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "Inferred Ct trajectory", 
       x = "Time (days since inferred exposure)", y = "Ct value")


p.predictive <- rbind(fit_dt_delta_pooled_summary, fit_dt_omicron_pooled_summary) %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = me, colour = voc)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = voc), alpha = 0.2) +
  scale_y_reverse() +
  custom_plot_theme(flip = TRUE) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "Inferred Ct trajectory", 
       x = "Time (days since inferred exposure)", y = "Ct value")  + 
  lims(y = c(45, 0))

# fit_dt_omicron_pooled_summary %>% 
#   ggplot() + 
#   geom_line(aes(x = as.numeric(time), y = me, colour = voc)) + 
#   # geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = voc), alpha = 0.2) +
#   scale_y_reverse() +
#   custom_plot_theme(flip = TRUE) + 
#   scale_fill_brewer(palette = "Set1") + 
#   labs(title = "Inferred Ct trajectory", 
#        x = "Time (days since inferred exposure)", y = "Ct value")

p.all <- p.all.posteriors/p.predictive
p.all
# ggsave("outputs/ct_trajectories.pdf",
#        p.all,
#        height = 6,
#        width = 12,
#        bg = "white")

