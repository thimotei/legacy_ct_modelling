dt_pop_priors_wide <- sample_pop_priors(5000,
                                        c_lod = 40,
                                        switch = FALSE, 
                                        data_format = "wide", 
                                        scale_type = "natural")

dt_pop_priors_long <- sample_pop_priors(1000000,
                                        c_lod = 40,
                                        switch = FALSE, 
                                        data_format = "long", 
                                        scale_type = "natural")

dt_pop_priors_long[!parameter %in% c("c_s", "t_s")] |> 
  ggplot() + 
  geom_density(aes(x = value, fill = parameter)) + 
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()


dt_pop_prior_sim <- simulate_cts(dt_pop_priors_wide,
                                 time_range = seq(0, 30, 0.1),
                                 obs_noise = FALSE, 
                                 by = c("t", "id"))

dt_pop_prior_sum <- dt_pop_prior_sim[,
                                     .(lo = quantile(ct_value, 0.025),
                                       me = quantile(ct_value, 0.5),
                                       hi = quantile(ct_value, 0.975)),
                                     by = c("t")]

dt_pop_prior_sim[exp_ct <= 40 & t > t_p | t <= t_p] |>
  ggplot() +
  geom_line(aes(x = t, y = exp_ct, group = id),
            alpha = 0.1) +
  scale_y_reverse() +
  theme_minimal()

dt_pop_prior_sum |>
  ggplot() +
  geom_line(aes(x = t, y = me)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = hi), alpha = 0.1) +
  scale_y_reverse() +
  theme_minimal()
