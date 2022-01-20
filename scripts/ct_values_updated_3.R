library(data.table)
library(lubridate)
library(patchwork)
library(rstan)
library(cowplot)
library(stringr)

source("R/stan_data.R")
source("R/custom_plot_theme.R")

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
  .[result == 0, ct_value := 45]

# dt.ct[is.na(t_symp_onset), t_symp_onset := t_symp_mean] %>% 
#   .[, t_first_test := t_first_test + 10, 
#                   by = c("id", "infection_id")] %>% 
#   .[, t_symp_onset := t_symp_onset + 10, 
#     by = c("id", "infection_id")] %>% 
#   .[, onset_rel := min(t_symp_onset), 
#     by = c("id", "infection_id")]

# dt.ct.delta[, onset_rel := as.numeric(symptom_onset - min(date), units = "days"), 
#                      by = c("id", "infection_id")] %>% 
#   # .[is.na(onset_rel) == TRUE, onset_rel := 3] %>% 
#   .[, onset_rel := onset_rel - min(onset_rel, na.rm = TRUE) + 10] %>% 
#   .[, t_first_test := as.numeric(date - min(date), units = "days"), 
#     by = c("id", "infection_id")] %>% 
#   .[, t_symp_onset := as.numeric(date - symptom_onset, units = "days"), 
#     by = c("id", "infection_id")] %>% 
#   .[, t_first_test := t_first_test -  min(onset_rel, na.rm = TRUE) + 10] %>% 
#   .[, t_symp_onset := t_symp_onset -  min(onset_rel, na.rm = TRUE) + 10]
#   

# dt.plot <- dt.ct[,  `:=`  (result = factor(result,
#                                            levels = c(1, 0), 
#                                            labels = c("Positive", "Negative")),
#                            voc = factor(voc,
#                                         levels = c("alpha", "delta", "omicron", "negative"),
#                                         labels = c("Alpha", "Delta", "Omicron", "Negative")))]


mn <- dt.ct[, min(ct_value, na.rm = TRUE)]
mx <- dt.ct[, max(ct_value, na.rm = TRUE)]

# scale Ct values
dt.ct[, ct_scaled := (ct_value - mn)/(mx - mn)]

# choosing with voc to fit to 
dt.ct.delta <- dt.ct[(voc == "omicron" | voc == "negative")] %>%
  .[, obs_total := .N, by = c("id", "infection_id")] %>%
  .[obs_total > 2] %>% 
  .[, onset_rel := as.numeric(symptom_onset - pmin(symptom_onset, date, na.rm = TRUE) + 10, units = "days"),
    by = id] %>% 
  .[, day_rel := as.numeric(date - pmin(symptom_onset, date, na.rm = TRUE) + 10, units = "days"),
    by = id] %>% 
  .[is.na(onset_rel) == FALSE] %>% 
  .[result == "Positive", result_num := 1] %>% 
  .[result == "Negative", result_num := 0] %>% 
  .[, id := .GRP, by = id]

stan_data_delta <- stan_data_fun(dt.ct.delta)

# dt.plot[(voc == "Delta" | voc == "Negative")] %>% 
#   .[, total_obs := .N, by = c("id", "infection_id")] %>% 
#   .[total_obs > 2] %>% 
#   ggplot() + 
#   geom_point(aes(x = t_first_test, y = ct_value,
#                  group = result,
#                  colour = result)) + 
#   geom_vline(aes(xintercept = onset_rel), linetype = "dashed") +
#   facet_wrap(~id, scales = "free")

# dt.ct[(voc == "delta" | voc == "negative") & is.na(symptom_onset) == FALSE] %>% 
#   .[, total_obs := .N, by = c("id", "infection_id")] %>% 
#   .[total_obs > 2] %>% 
#   ggplot() + 
#   geom_point(aes(x = t_symp_onset, y = ct_value,
#                  group = interaction(id, result), colour = result)) + 
#   facet_wrap(~id, scales = "free")


options(mc.cores = parallel::detectCores()) 
mod <- rstan::stan_model("stan/ct_trajectory_model.stan")

fit_delta <- sampling(mod,
                      chains = 4,
                      data = stan_data_delta,
                      iter = 2000)

res <- extract(fit_delta)

dt.tmp <- res$c_p_mean %>% 
  melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn]

omicron_c_peak <- dt.tmp %>% 
  ggplot() +
  geom_density(aes(smp), fill = "lightblue") + 
  xlim(0, 45) +
  labs(x = "Ct value", y = "Probability",
       title = "Posterior distribution for peak Ct value") +
  custom_plot_theme()

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

tmp <- data.table::melt(as.data.table(res$c_p)[, iterations := 1:.N],
                        id.vars = c("iterations"),
                        value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
                        ][, num_id := as.numeric(num_id)] %>%
  .[, smp := (mx - mn) * smp + mn]

tmp %>% 
  ggplot() + 
  geom_density(aes(x = smp))

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

tmp
