library(rstan)
library(stringr)

source("R/stan_data_simulated.R")

stan_data_simulated <- stan_data_simulated(100)

mod_sim <- stan_model("stan/ct_trajectory_simulate_data.stan")
fit_sim<- sampling(mod_sim,
                   chains = 100,
                   data = stan_data_simulated,
                   iter = 1,
                   algorithm = "Fixed_param")

res_sim <- extract(fit_sim) 

ct.dt <- res_sim$ct %>% 
  data.table() %>% 
  melt(., measure.vars = patterns("V"),
       variable.name = "time",
       value.name = "ct") %>% 
  .[, id := 1:100, by = "time"] %>% 
  .[, time := as.numeric(str_remove(time, "V"))] %>% 
  .[, ct := (mx - mn) * ct + mn] %>% 
  .[order(id, time)]

ct.dt[time < 30 & id %in% 1:20] %>% 
  ggplot() + 
  geom_line(aes(x = time, y = ct, group = id, colour = factor(id))) + 
  theme(legend.position = "none")

