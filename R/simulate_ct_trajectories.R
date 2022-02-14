simulate_ct_trajectories <- function(t_max, t_stepsize,
                                     cp_min, cp_max, 
                                     cs_min, cs_max, 
                                     te_min, te_max, 
                                     tp_min, tp_max, 
                                     ts_min, ts_max, 
                                     tlod_min, tlod_max,
                                     sigma_obs) {
  out_dt <- data.table(id = sort(rep(1:n, t_max/t_stepsize)), 
                       t = rep(t_input/t_stepsize, n), 
                       c0 = sort(rep(rep(c0, n), t_max/t_stepsize)),
                       cp = sort(rep(runif(n, cp_min, cp_max), t_max/t_stepsize)),
                       cs = sort(rep(runif(n, cs_min, cs_max), t_max/t_stepsize)),
                       te = sort(rep(runif(n, te_min, te_max), t_max/t_stepsize)),
                       tp = sort(rep(runif(n, tp_min, tp_max), t_max/t_stepsize)),
                       ts = sort(rep(runif(n, ts_min, ts_max), t_max/t_stepsize)),
                       tlod = sort(rep(runif(n, tlod_min, tlod_max), t_max/t_stepsize))) %>% 
    .[, ct_value := ct_hinge_long(t, c0 = c0, cp = cp, cs = cs,
                                clod = clod, te = te, tp = tp, ts = ts, 
                                tlod = tlod),
      by = c("id", "t")] %>% 
    .[, ct_value_noisey := rtruncnorm(1, b = 40, mean = ct_value, sd = sigma_obs),
      by = c("id")] %>%
    .[, t := as.numeric(t)] %>% 
    .[(t < te + 1 | t > te + tp + ts + tlod + 1), pcr_res := 0] %>%
    .[t > te & t < te + tp + ts + tlod + 1, pcr_res := 1] %>%
    .[, pcr_res := factor(pcr_res)] %>% 
    .[, ct_value_std := (ct_value_noisey - min(ct_value_noisey))/(max(ct_value_noisey) - min(ct_value_noisey))]
  
  return(out_dt)
}
