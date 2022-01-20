library(data.table)
library(ggplot2)
library(rstan)

dt.ct <- fread("data/ct_values_dry_vs_wet.csv") %>% 
  setnames(.,
           c("Study number", "sample date", "barcode",
                "swab type", "Ct", "VOC"),
           c("id", "date", "barcode", "swab_type", "ct_value",
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
  .[, ct_value := as.numeric(ct_value)] %>% 
  .[is.na(ct_value), ct_value := NA] %>% 
  .[result == 0, ct_value := 45]


t1 <- dt.ct %>% 
  .[, obs := .N, by = c("id", "date")] %>%  
  .[obs > 1] %>% 
  .[, c("id", "date", "swab_type", "ct_value")]


t2 <- t1[, .SD[swab_type != shift(swab_type)], by = id] %>% 
  .[order(id, date)] %>% 
  .[, obs:= .N, by = c("id", "date")] %>% 
  .[obs > 1] %>% 
  na.omit()

t3 <- t2[swab_type == "VTM"]
t4 <- t2[swab_type == "Dry"]
t5 <- t2[swab_type == "pipeline"][, swab_type := "Dry"]

t6 <- rbind(t4, t5) 

t7 <- merge.data.table(t3, t6, by = c("id", "date"))

#--- frequentist linear regression fits

p.vtm.vs.wet.frequentist <- t7[ct_value.x < 45] %>% 
  ggplot(aes(x = ct_value.x, y = ct_value.y)) + 
  geom_point() +
  labs(x = "Ct value of VTM test", y = "Ct value of Dry swab",
       title = "Frequentist fit") + 
  geom_smooth(method = "lm") +
  theme(legend.position = "none") +
  ylim(13, 30)

# ggsave("outputs/dry_vs_wet_scatter_plot.png",
#        p1,
#        width = 6,
#        height = 6,
#        bg = "white")

cor.test(t7[ct_value.x < 45]$ct_value.x, t7[ct_value.x < 45]$ct_value.y, method = "pearson")

#--- Bayesian linear regression fits

stan_data <- list(
  N = t7[ct_value.x < 45, .N],
  x = t7[ct_value.x < 45, ct_value.x],
  y = t7[ct_value.x < 45, ct_value.y])

options(mc.cores = parallel::detectCores())
mod <- stan_model("stan/dry_vs_wet_linear_regression.stan")
fit <- sampling(mod, 
                data = stan_data,
                chains = 4,
                iter = 2000,
                warmup = 1000,
                control = list(adapt_delta = 0.99))   

res <- extract(fit)

# a data.table summarising the sampled values into a median and 95% credible
# intervals
y_pred_dt <- melt(data.table(res$y_pred),
                  measure = patterns("V"),
                  variable.name = "iteration") %>%
  .[, iteration := stringr::str_remove(string = iteration, pattern = "V")] %>%
  .[, iteration := as.numeric(iteration)] %>%
  .[, .(median = median(value),
        UQ = quantile(value, 0.975),
        LQ = quantile(value, 0.025)), by = "iteration"] %>% 
  .[, x_data:= stan_data$x] %>% 
  .[, y_data := stan_data$y] %>% 
  setcolorder(., c("iteration", "x_data", "y_data"))

# plotting the median and 95% CrI of the fitted curve and the empirical
# (simulated) data together
p.vtm.vs.wet.bayesian <- y_pred_dt %>% 
  ggplot(aes(x = x_data)) +
  geom_line(aes(y = median)) + 
  geom_ribbon(aes(ymin = LQ, ymax = UQ), fill = "dodgerblue", alpha = 0.2) + 
  geom_point(aes(x = x_data, y = y_data), inherit.aes = FALSE) +
  #geom_smooth(method = "lm", aes(x = x_data, y = y_data)) +
  scale_fill_identity(name = "95% confidence/credible intervals",
                      guide = "legend", labels = c("Bayesian fit", "Frequentist fit")) +
  labs(x = "Ct value of VTM test", y = "Ct value of Dry swab",
       title = "Bayesian fit")

p.together <- p.vtm.vs.wet.frequentist + p.vtm.vs.wet.bayesian

# ggsave("outputs/vtm_vs_dry.png",
#        p.together,
#        width = 8, 
#        height = 4,
#        bg = "white")

# investigating the posterior distributions of the parameters of the curve
# we can also look at the posterior predictive distributions, by looking at
# y_pred, but there are 100 of them and we can equally investigate their
# goodness of fit by looking at the fitted curve
stan_dens(fit, pars = c("alpha", "beta", "sigma"))

beta.dt <- data.table(reshape2::melt(res$beta, 
                                     measure = patterns("V"),
                                     variable.name = "iteration"))


adjustment.mi <- quantile(beta.dt[, value], 0.5)
adjustment.lo <- quantile(beta.dt[, value], 0.025)
adjustment.hi <- quantile(beta.dt[, value], 0.975)

adjustment.dt <- data.table(mi = adjustment.mi, 
                            lo = adjustment.lo,
                            hi = adjustment.hi)
