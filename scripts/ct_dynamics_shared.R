library(data.table)
library(ggplot2)
library(lubridate)
library(cowplot)
library(patchwork)
library(encryptr)

#--- defining a custom plot theme for ggplot objects
custom_plot_theme <- function(flip = FALSE, legend_arg = FALSE) {
  
  custom_plot_theme = list(
    theme_cowplot(11),
    theme(strip.placement = "outside",
          strip.background = element_blank()),
    geom_vline(aes(xintercept = -Inf)),
    guides(color = guide_legend(override.aes = list(fill = NA))))
  
  if(flip == FALSE) { 
    append(custom_plot_theme, geom_hline(aes(yintercept = -Inf)))
  }
  
  if(legend_arg == FALSE) {
    append(custom_plot_theme,
           theme(legend.title = element_blank()))
  }
  
  return(custom_plot_theme)
}

#--- defining function for dry vs wet swab adjustment 
#--- we adjust the VTM swabs downwards
adjustment.fun <- function(alpha, beta, x) {alpha + beta*x}
adjustment.dt <- data.table(param = c("alpha", "beta"),
                            me = c(0.949591208744795, 0.807592348470527),
                            lo = c(-1.857339823018790, 0.685903357311594),
                            hi = c(4.02657501060788, 0.919626934297874))

#--- loading in encrypted R object
load("~/Dropbox/Legacy_LSHTM/ct_dynamics.Rdata") 

dt.ct <- dt.ct %>% 
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
  .[swab_type == "VTM" & result == "Negative", ct_adjusted := 40] %>%
  .[swab_type == "VTM" & result != "Negative", ct_adjusted := adjustment.fun(adjustment.dt[param == "alpha", me],
                                                                             adjustment.dt[param == "beta", me],
                                                                                           ct)] %>%
  .[swab_type != "VTM", ct_adjusted := ct]

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

#--- counting number of positive swabs by individual
no_pos_cts <- dt.ct[(result == "Positive" | result == "Inconclusive") & is.na(ct_adjusted) == FALSE] %>% 
  .[, .N, by = c("ID", "infection_ID")]

dt.ct <- merge.data.table(dt.ct, no_pos_cts, by = c("ID", "infection_ID")) %>%  
  setnames(., 
           c("N", "number_vaccines (14 days pre ix)"), 
           c("no_pos_results", "no_vaccines"))

#--- Plotting results, first row - "with negative results"
p1.all <- dt.ct %>% 
  .[days_since_first_test < 20 & (VOC == "Delta" | VOC == "Omicron") & result != "Invalid" & result != ""] %>%
  .[, .SD[any(result == "Negative")], by = ID] %>%
  ggplot(aes(x = days_since_first_test, y = ct, colour = interaction(factor(VOC)))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme() + 
  labs(x = "Time (days)",
       y = "Ct value",
       title = "With negative test results",
       subtitle = "Since first test",
       colour = "VOC") + 
  scale_colour_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  lims(x = c(-5, 20), y = c(42, 10))

p2.all <- dt.ct %>% 
  .[days_since_symptom_onset > -4 & days_since_symptom_onset < 20 & (VOC == "Delta" | VOC == "Omicron") & result != "Invalid" & result != ""] %>% 
  ggplot(aes(x = days_since_symptom_onset, y = ct, colour = factor(VOC))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme() + 
  labs(x = "Time (days)",
       y = "Ct value",
       subtitle = "Since symptom onset",
       colour = "VOC") + 
  scale_colour_brewer(palette = "Set1") + 
  lims(x = c(-5, 20), y = c(42, 10)) + 
  theme(plot.subtitle = element_text(hjust = 0.5))

p.all <- p1.all + p2.all

#--- Plotting results, second row - "without negative results"
p1.no.neg <- dt.ct %>% 
  .[days_since_first_test < 20 & (VOC == "Delta" | VOC == "Omicron") & (result == "Positive" | result == "Inconclusive")] %>% 
  ggplot(aes(x = days_since_first_test, y = ct, colour = factor(VOC))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme() + 
  labs(x = "Time (days)",
       y = "Ct value",
       title = "Without negative test results",
       colour = "VOC") + 
  scale_colour_brewer(palette = "Set1") +
  lims(y = c(42, 10)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  lims(x = c(-5, 20))

p2.no.neg <- dt.ct %>% 
  .[days_since_symptom_onset > -4 & days_since_symptom_onset < 20 & (VOC == "Delta" | VOC == "Omicron") & (result == "Positive" | result == "Inconclusive")] %>% 
  ggplot(aes(x = days_since_symptom_onset, y = ct, colour = factor(VOC))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme() + 
  labs(x = "Time (days)", y = "Ct value", colour = "VOC") + 
  scale_colour_brewer(palette = "Set1") + 
  lims(x = c(-5, 20), y = c(42, 10)) +
  theme(plot.subtitle = element_text(hjust = 0.5))

p.all.no.neg <- p1.no.neg + p2.no.neg 
p.all.both <- p.all/p.all.no.neg + plot_layout(guides = "collect") 

#--- Plotting results, third row - "stratifying omicron curves by numbers of vaccine"
p1.vaccines <- dt.ct %>% 
  .[days_since_first_test < 20 & (VOC == "Omicron") & result != "Invalid" & result != ""] %>% 
  .[, .SD[any(result == "Negative")], by = ID] %>%
  ggplot(aes(x = days_since_first_test, y = ct, colour = interaction(factor(no_vaccines)))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme(legend_arg = TRUE) + 
  labs(x = "Time (days)", y = "Ct value",
       title = "Omicron curves stratified by number of vaccines",
       subtitle = "Since first test",
       colour = "Number of \n vaccines") + 
  scale_colour_brewer(palette = "Set2") +
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  lims(x = c(-5, 20)) + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

p2.vaccines <- dt.ct %>% 
  .[days_since_symptom_onset > -4 & days_since_symptom_onset < 20 & (VOC == "Omicron") & result != "Invalid" & result != ""] %>% 
  ggplot(aes(x = days_since_symptom_onset, y = ct, colour = factor(no_vaccines))) + 
  geom_jitter(width = 0, height = 0.5, alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  stat_smooth(alpha = 0.15) + 
  scale_y_reverse() + 
  custom_plot_theme(legend_arg = TRUE) + 
  labs(x = "Time (days)",
       y = "Ct value",
       subtitle = "Since symptom onset",
       colour = "Number of \n vaccines") + 
  scale_colour_brewer(palette = "Set2") + 
  lims(y = c(42, 10)) + 
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  lims(x = c(-5, 20)) + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

p.vaccines <- p1.vaccines + p2.vaccines + plot_layout(guides = "collect")
p.all.vaccines <- p.all.both/p.vaccines

#--- saving overall plot with all three rows
# ggsave("outputs/ct_dynamics_adjusted_vaccines.png",
#        p.all.vaccines,
#        width = 7,
#        height = 9,
#        bg = "white")

