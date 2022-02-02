library(data.table)
library(lubridate)
library(patchwork)
library(rstan)
library(cowplot)
library(stringr)
library(purrr)

source("R/custom_plot_theme.R")

dt.raw <- fread("data/ct_values_clean.csv") 

dt.ct <- fread("data/ct_values_clean.csv") %>% 
  .[, swab_date := dmy(swab_date)] %>% 
  .[barcode %like% "49U", swab_date := swab_date - 1] %>% 
  .[, ct := as.numeric(ct)]

dt.reduced <- dt.ct[,  c("ID", "barcode", "swab_date", "swab_type",
                         "ct")]

dt.vtm <- dt.reduced[swab_type == "VTM"]
dt.dry <- dt.reduced[swab_type == "Dry"]

dt.both <- merge.data.table(dt.vtm, dt.dry, by = c("ID", "swab_date")) %>% 
  setnames(., c("ct.x", "ct.y"), c("ct.vtm", "ct.dry"))

stan_data <- list(
  N = dt.both[ct.vtm < 45 & is.na(ct.vtm) == FALSE & is.na(ct.dry) == FALSE, .N],
  x = dt.both[ct.vtm < 45 & is.na(ct.vtm) == FALSE & is.na(ct.dry) == FALSE, ct.vtm],
  y = dt.both[ct.dry < 45 & is.na(ct.vtm) == FALSE & is.na(ct.dry) == FALSE, ct.dry])

options(mc.cores = parallel::detectCores())
mod <- stan_model("stan/dry_vs_wet_linear_regression.stan")
fit <- sampling(mod, 
                data = stan_data,
                chains = 4,
                iter = 2000,
                warmup = 1000,
                control = list(adapt_delta = 0.99))   

res <- extract(fit)

p.vtm.vs.wet.frequentist <- dt.both[ct.vtm < 45] %>% 
  ggplot(aes(x = ct.vtm, y = ct.dry)) + 
  geom_point() +
  labs(x = "Ct value of VTM test", y = "Ct value of Dry swab",
       title = "Frequentist fit") + 
  geom_smooth(method = "lm") +
  theme(legend.position = "none")

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
# stan_dens(fit, pars = c("alpha", "beta", "sigma"))
alpha.dt <- data.table(reshape2::melt(res$alpha, 
                                      measure = patterns("V"),
                                      variable.name = "iteration")) %>% 
  .[, .(param = "alpha", 
        me = quantile(value, 0.5), 
        lo = quantile(value, 0.025),
        hi = quantile(value, 0.975))]

beta.dt <- data.table(reshape2::melt(res$beta, 
                                     measure = patterns("V"),
                                     variable.name = "iteration")) %>% 
  .[, .(param = "beta",
        me = quantile(value, 0.5), 
        lo = quantile(value, 0.025),
        hi = quantile(value, 0.975))]

param.dt <- rbind(alpha.dt, beta.dt)

fwrite(param.dt, "data/adjustment_params.csv")
