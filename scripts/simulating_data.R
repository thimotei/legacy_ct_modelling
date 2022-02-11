library(data.table)
library(ggplot2)
library(rstan)
library(stringr)
library(purrr)

source("R/stan_data_simulated.R")
source("R/get_inc_period.R")

ct.data <- fread("data/ct_values_clean.csv")

mx <- ct.data %>% 
  .[ct != "" & ct != "unknown", max(ct, na.rm = TRUE) %>% as.numeric]

mn <- ct.data %>% 
  .[ct != "", min(ct, na.rm = TRUE) %>% as.numeric]

stan_data_simulated <- stan_data_simulated(P_arg = 100, t_max = 40)

# replicating ordered parameters in Stan, making sure all timing values are well-ordered
init_fun <- function(n.chains) rep(list(
  t_e = runif(1, 1, 5),
  t_p = runif(1, 5, 7),
  t_s = runif(1, 7, 10),
  t_lod = runif(1, 10, 13)),
  n.chains)

# mod_sim <- stan_model("stan/ct_trajectory_simulate_data.stan")
mod_sim <- stan_model("stan/ct_trajectory_simulate_data_non_hierarchical.stan")
fit_sim <- sampling(mod_sim,
                   chains = 100,
                   data = stan_data_simulated,
                   warmup = 0,
                   iter = 1,
                   init = init_fun,
                   algorithm = "Fixed_param")

res_sim <- extract(fit_sim) 

ct.dt <- res_sim$ct %>% 
  data.table() %>% 
  melt(., measure.vars = patterns("V"),
       variable.name = "time",
       value.name = "ct") %>% 
  .[, id := 1:100, by = "time"] %>% 
  .[, time := as.numeric(str_remove(time, "V"))] %>% 
  .[, ct_original := (mx - mn) * ct + mn] %>% 
  .[order(id, time)] 

mx_adj <- ct.dt[, max(ct)]
mn_adj <- ct.dt[, min(ct)]

ct.dt.samp <- ct.dt %>% 
  .[, .SD[time %in% sample(.N, sample(3:6, 1))], by = "id"] %>% 
  setnames(., "time", "t_test") %>% 
  .[order(id, t_test)] %>% 
  .[ct == mx_adj, pcr_res := 0] %>% 
  .[ct < mx_adj, pcr_res := 1] %>% 
  .[ct < mx_adj, ct_samp := as.numeric(truncnorm::rtruncnorm(1, a = mn_adj, b = mx_adj, mean = ct, sd = 0.1))] %>% 
  .[ct == mx_adj, ct_samp := ct] 

# choosing subset to fit to
ct.dt.samp.stan <- ct.dt.samp[id %in% 1:20] 

# plotting subset
ct.dt.samp.stan %>%
  ggplot() +
  geom_line(aes(x = t_test, y = ct_samp)) +
  geom_point(aes(x = t_test, y = ct_samp)) +
  facet_wrap(~id)

#--- attempting to fit simulated data

stan_data <- list(
  P = ct.dt.samp.stan[, uniqueN(id)],
  N = ct.dt.samp.stan[, .N],
  ct_value = ct.dt.samp.stan[, ct],
  id = ct.dt.samp.stan[, id],
  day_rel = ct.dt.samp.stan[, t_test],
  pcr_res = ct.dt.samp.stan[, pcr_res],
  t_e = 0,
  c_0 = (40 - mn)/(mx - mn),
  c_lod = (40 - mn)/(mx - mn),
  lmean = get_inc_period()$inc_mean_p[1],
  lsd = get_inc_period()$inc_sd_p[2]
)

mod <- stan_model("stan/ct_trajectory_model_individual.stan")
options(mc.cores = parallel::detectCores())
fit <- sampling(mod, data = stan_data)

#--- inspecting fits
res <- extract(fit)

dt.t.e <- res$T_e %>% 
  reshape2::melt() %>% 
  data.table()  %>% 
  setnames(., "Var2", "id")

dt.t.peak <- res$t_p_mean %>% 
  reshape2::melt() %>% 
  data.table()

dt.t.switch <- res$t_s_mean %>% 
  reshape2::melt() %>% 
  data.table()

dt.t.lod <- res$t_lod_mean %>% 
  reshape2::melt() %>% 
  data.table()

dt.c.peak <- res$c_p_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn]

dt.c.switch <- res$c_s_mean %>% 
  reshape2::melt() %>% 
  data.table() %>% 
  .[, smp := (mx - mn) * value + mn]

T.e.posteriors <- dt.t.e %>%
  ggplot() +
  geom_density(aes(x = value), alpha = 0.3) +
  facet_wrap(~id, scales = "free") + 
  xlim(0, 10) +
  labs(x = "Time", y = "Probability",
       title = "Posterior distribution for exposure time") +
  custom_plot_theme() +
  scale_fill_brewer(palette = "Set2") 

dt.c.p.s <- rbind(dt.c.peak[, type := "peak"],
                  dt.c.switch[, type := "switch"])

dt.c.p.s %>% 
  ggplot() + 
  geom_density(aes(x = smp, fill = type), alpha = 0.4) + 
  xlim(0, 40)

dt.t.p.s.lod <- rbind(dt.t.peak[, type := "peak"],
                  dt.t.switch[, type := "switch"],
                  dt.t.lod[, type := "lod"])

dt.t.p.s.lod %>% 
  ggplot() + 
  geom_density(aes(x = value, fill = type), alpha = 0.4) + 
  xlim(0, 10)


# ggsave("outputs/all_te_posteriors.png", 
#        T.e.posteriors,
#        width = 10,
#        height = 10,
#        bg = "white")


ct_melt_fun <- function(id) {
  
  dt.out <- res$ct[, id,] %>%
    reshape2::melt() %>% 
    data.table()
  
  return(dt.out)
  
}

dt.tmp <- 1:stan_data$P %>% 
  map_df(~ct_melt_fun(.)) %>% 
  setnames(., 
           c("Var1", "Var2", "value"),
           c("iteration", "time", "ct"))

dt.predictive <- 1:stan_data$P %>% 
  map_df(~ct_melt_fun(.)) %>% 
  setnames(., 
           c("Var1", "Var2", "value"),
           c("iteration", "time", "ct")) %>% 
  .[, ct := (mx - mn) * ct + mn] %>% 
  .[order(time)] %>%  
  .[, id := 1:stan_data$P,
    by = c("iteration", "time")] %>% 
  .[, .(me = quantile(ct, 0.5), 
        lo = quantile(ct, 0.025),
        hi = quantile(ct, 0.975)),
    by = c("id", "time")]

dt.predictive[id %in% 1:10] %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = me)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_point(data = ct.dt.samp[id %in% 1:10], aes(x = t_test, y = ct_original)) +
  facet_wrap(~id)

dt.predictive[id %in% 1:10] %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = me)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) 

p.all <- p.all.posteriors/p.predictive
p.all
# ggsave("outputs/ct_trajectories.pdf",
#        p.all,
#        height = 6,
#        width = 12,
#        bg = "white")



