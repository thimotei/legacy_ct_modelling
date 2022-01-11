library(data.table)
library(lubridate)
library(rstan)
library(patchwork)
library(cowplot)

source("R/ct_hinge_funs.R")
source("R/stan_data.R")

dt.raw <- fread("data/pcr_data/serial_ct_infection_history.csv")

setnames(dt.raw, 
         c("Study number", "sample date", "barcode",
           "Ct", "sample order", "Day post swab 1", "VOC"), 
         c("id", "date", "sample_barcode", "ct_value", "obs_counter",
           "day_rel", "voc"))

dt.raw[, date := dmy(date)]
#dt.raw <- dt.raw[!(id == 977 & date == "2021/11/15")]
#dt.raw <- dt.raw[is.na(day_rel) == FALSE]
dt.raw[voc == "Delta", voc := "delta"]
dt.raw[voc == "Omicron", voc := "omicron"]
dt.raw[voc == "OMicron", voc := "omicron"]
dt.raw[voc != "negative", result := 1]
dt.raw[voc == "negative", `:=` (result = 0, voc = NA)]
dt.raw[is.na(ct_value) == TRUE, ct_value := 45]
#dt.raw[result == 1, day_rel := date = min(date), by = id]
dt.raw[, obs_counter := NULL]
dt.raw[, obs_both := 1:.N, by = id]
dt.raw[, obs_max := max(obs_both), by = id]
dt.raw[result == 1, obs_pos := 1:.N, by = id]
dt.raw[result == 0, obs_neg := 1:.N, by = id] 

dt.raw[, infection_id := .GRP, by = c("id", "infection_number")]
dt.raw[is.na(voc) == TRUE, voc := "negative"]

# leaving out participant 969 for the time being, as the timeline of
# PCR results are dubious
dt.raw <- dt.raw[id != 969]
dt.raw <- dt.raw[voc != "alpha"]

# p.all.raw <- dt.raw[, .SD[.N > 1], by = infection_id] %>% 
#   ggplot() + 
#   geom_point(aes(x = date, y = ct_value, colour = factor(voc))) + 
#   geom_line(aes(x = date, y = ct_value)) + 
#   facet_wrap(~infection_id, scales = "free") +
#   scale_y_reverse() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   scale_x_date(date_labels = "%d-%b")
# 
# ggsave("all_raw_data.png",
#        p.all.raw,
#        width = 12,
#        height = 12,
#        bg = "white")
# 
# p.delta.raw <- dt.raw[, .SD[.N > 1], by = infection_id] %>% 
#   .[voc == "delta" | voc == "negative"] %>% 
#   .[, .SD[any(voc == "delta")], by = infection_id] %>% 
#   ggplot() + 
#   geom_point(aes(x = date, y = ct_value, colour = factor(voc))) + 
#   geom_line(aes(x = date, y = ct_value)) + 
#   facet_wrap(~infection_id, scales = "free") +
#   scale_y_reverse() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   scale_x_date(date_labels = "%d-%b")
# 
# ggsave("delta_raw_data.png",
#        p.delta.raw,
#        width = 8,
#        height = 8,
#        bg = "white")
# 
# p.omicron.raw <- dt.raw[, .SD[.N > 1], by = infection_id] %>% 
#   .[voc == "omicron" | voc == "negative"] %>% 
#   .[, .SD[any(voc == "omicron")], by = infection_id] %>% 
#   ggplot() + 
#   geom_point(aes(x = date, y = ct_value, colour = factor(voc))) + 
#   geom_line(aes(x = date, y = ct_value)) + 
#   facet_wrap(~infection_id, scales = "free") +
#   scale_y_reverse() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   scale_x_date(date_labels = "%d-%b")
# 
# ggsave("omicron_raw_data.png",
#        p.omicron.raw,
#        width = 8,
#        height = 8,
#        bg = "white")


# adding artificial negative results for 5 days prior to the first swab recorded
# to help inference process

dt.symptom <- fread("data/pcr_data/symptom_onset.csv") %>% 
  setnames(.,
         c("participant_id", "symptom_onset"),
         c("id", "date_onset")) %>% 
  .[, date_onset := dmy(date_onset)] %>% 
  .[is.na(date_onset) == FALSE]

dt.all <- merge.data.table(dt.raw, dt.symptom, all.x = TRUE,
                           by = c("id", "infection_number"))


p1 <- dt.raw[, .SD[.N > 1], by = infection_id] %>%
  .[voc == "delta" | voc == "negative"] %>%
  .[, .SD[any(voc == "delta")], by = infection_id] %>%
  ggplot() +
  geom_point(aes(x = day_rel, y = ct_value, colour = factor(voc))) +
  #geom_line(aes(x = day_rel, y = ct_value)) +
  geom_smooth(aes(x = day_rel, y = ct_value)) +
  scale_y_reverse() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_cowplot(11) + 
  xlab("Time since first test") + 
  ylab("Ct value") +
  scale_color_manual(labels = c("Positive", "Negative"), 
                     values = c("blue", "red")) +
  labs(title = "Delta")

p2 <- dt.raw[, .SD[.N > 1], by = infection_id] %>%
  .[voc == "omicron" | voc == "negative"] %>%
  .[, .SD[any(voc == "omicron")], by = infection_id] %>%
  ggplot() +
  geom_point(aes(x = day_rel, y = ct_value, colour = factor(voc))) +
  #geom_line(aes(x = day_rel, y = ct_value)) +
  geom_smooth(aes(x = day_rel, y = ct_value)) +
  scale_y_reverse() +
  theme_cowplot(11) + 
  xlab("Time since first test") + 
  ylab("Ct value") +
  theme(axis.text.x = element_text( hjust = 1), legend.position = "none") +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Omicron")

p1 + p2 + plot_layout(guides = 'collect')

dt.all[result == 1, first_pos_date := min(date), by = c("id", "infection_number")]

dt.all %>% 
  ggplot(aes(x = first_pos_date, y = date_onset)) + 
  geom_point() 

# preparing stan data
mn <- dt.raw[, min(ct_value, na.rm = TRUE)]
mx <- dt.raw[, max(ct_value, na.rm = TRUE)]

# scale Ct values
dt.all[, ct_scaled := (ct_value - mn)/(mx - mn)]

dt.all[voc == "delta" | voc == "omicron" | voc == "alpha", result := 1]
dt.all[voc == "negative", result := 0]
dt.all <- dt.all[voc != "retest/inconsistent"]

# choosing with voc to fit to 
dt.together.both <- dt.all[, .SD[.N > 1], by = infection_id] %>% 
  .[voc == "omicron" | voc == "delta" | voc == "negative"] %>% 
  .[, .SD[any(voc == "omicron") | any(voc == "delta")], by = infection_id] %>% 
  .[, id := .GRP, by = id] %>% 
  .[, date_onset_rel := as.numeric(date_onset - min(date), units = "days"), by = id]

dt.together.delta <- dt.all[, .SD[.N > 1], by = infection_id] %>% 
  .[voc == "delta" | voc == "negative"] %>% 
  .[, .SD[any(voc == "delta")], by = infection_id] %>% 
  .[, id := .GRP, by = id] %>% 
  .[, date_onset_rel := as.numeric(date_onset - min(date), units = "days"), by = id]

dt.together.omicron <- dt.all[, .SD[.N > 1], by = infection_id] %>% 
  .[voc == "omicron" | voc == "negative"] %>% 
  .[, .SD[any(voc == "omicron")], by = infection_id] %>%
  .[, id := .GRP, by = id] %>% 
  .[, date_onset_rel := as.numeric(date_onset - min(date), units = "days"), by = id]

# Construct stan data
stan_data_delta <- stan_data_fun(dt.together.delta)
stan_data_omicron <- stan_data_fun(dt.together.omicron)
stan_data_both <- stan_data_fun(dt.together.both)

options(mc.cores = parallel::detectCores())
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

res_omicron <- extract(fit_omicron)

fit_omicron_dt <- as.data.frame(fit_omicron, pars = "ct") %>%
  data.table() %>%
  melt() %>% 
  .[, c("id", "time") := tstrsplit(variable, ",")] %>% 
  .[, id := str_remove(id, "ct\\[")] %>% 
  .[, time := str_remove(time, "]")] %>%
  .[, iteration := 1:.N, by = "variable"] %>% 
  .[, variable := NULL] %>% 
  .[, time := as.numeric(time)]

fit_omicron_dt_summary <- fit_omicron_dt[, .(me = quantile(value, 0.5), 
                                             lo = quantile(value, 0.025),
                                             hi = quantile(value, 0.975)),
                                         by = c("id", "time")]

tmp <- fit_omicron_dt_summary %>% 
  ggplot() + 
  geom_line(aes(x = as.numeric(time), y = lo)) + 
  geom_ribbon(aes(x = time, ymin = lo, ymax = hi, fill = "dodgerblue"), alpha = 0.1) +
  facet_wrap(~id)

ggsave("all_ct_fits.png",
       tmp,
       width = 10,
       height = 10,
       bg = "white")

df <- data.frame(c("5", "3", "8 [3 - 5]")) 
df2 = str_split_fixed(string = df[,1], pattern = "\\[", n = 2)               
df2[,2] = gsub(pattern = "\\]", replacement = "", x = df2[,2])

t_e = 0
c_0 = (45 - mn)/(mx - mn)
c_lod = (45 - mn)/(mx - mn)

big_out <- data.table(t_peak = exp(res_omicron$t_p_mean),
                      t_switch = exp(res_omicron$t_s_mean),
                      t_lod = exp(res_omicron$t_lod_mean),
                      c_switch = stan_data_omicron$c_lod * boot::inv.logit(res_omicron$c_s_mean),
                      c_peak = stan_data_omicron$c_lod * boot::inv.logit(res_omicron$c_s_mean) * boot::inv.logit(res_omicron$c_p_mean))


big_list <- list()
for(i in 1:nrow(big_out)){
  print(i)
  big_list[[i]] <- ct_hinge_dt(c_peak = big_out$c_p[i],
                               t_peak = big_out$t_p[i],
                               c_switch = big_out$c_s[i],
                               t_switch = big_out$t_s[i],
                               t_LOD = big_out$t_lod[i] + big_out$t_s[i] + big_out$t_p[i])
  
}

big.out.dt <- big_out[, ct := ct_hinge_dt(c_peak, t_peak, c_switch, t_switch, t_lod + t_switch + t_peak)]

ct_hinge_dt(big_out$c_peak, big_out$t_peak, big_out$c_switch, big_out$t_switch,
            big_out$t_lod + big_out$t_switch + big_out$t_peak)


vapply(X = seq(0, 40, 0.1), FUN = ct_hinge, FUN.VALUE = 1, c_zero = c_zero,
       c_switch = c_switch, t_switch = t_switch, t_eclipse = 0,
       c_LOD = c_lod, c_peak = c_peak, t_peak = t_peak, t_LOD = t_lod)

ct_hinge <- function(a, c_zero, 
                     c_peak, c_switch, 
                     c_LOD, t_eclipse, 
                     t_peak, t_switch, t_LOD)

res_omicron$ct[,,1] %>% data.table() %>% melt(variable.name = "id") %>% 
  .[, median(value), by = id]


dim1 <- c('iter')
dim2 <- c('iter', 'id')
dim3 <- c('iter', 'id', 'time')

cbind(expand.grid(dim1,  dim2, dim3), val = as.vector(data))

cbind(expand.grid(dim1, dim2, dim3), val = as.vector(res_omicron$ct))

expand.grid(res_omicron)


melt(data.table(res_omicron))


