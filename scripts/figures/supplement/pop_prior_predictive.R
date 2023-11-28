source("scripts/setup/main.R")

# sampling from the covariate-level priors
dt_pop_priors_wide <- sample_pop_priors(100000,
                                        c_lod = 40,
                                        switch = FALSE, 
                                        data_format = "wide", 
                                        scale_type = "transformed")

# transforming the parameters to the R code model scale
dt_pop_priors_wide_mod <- transform_to_model(dt_pop_priors_wide)

# simulating the trajectories
dt_pop_prior_sim <- simulate_cts(
  dt_pop_priors_wide_mod,
  time_range = seq(0, 30, 0.5),
  obs_noise = FALSE, 
  by = c("t", "id"))

# summarising the trajectories
dt_pop_prior_sum <- dt_pop_prior_sim[
  , .(lo = quantile(ct_value, 0.025),
      me = quantile(ct_value, 0.5),
      hi = quantile(ct_value, 0.975)),
  by = c("t")]


# plotting the median and 95% credible intervals
p_prior_predictive <- dt_pop_prior_sum |>
  ggplot() +
  geom_line(aes(x = t, y = me), linetype = "dashed") +
  geom_ribbon(aes(x = t, ymin = lo, ymax = hi), 
              fill = "dodgerblue",
              colour = "dodgerblue",
              alpha = 0.1) +
  scale_y_reverse() +
  labs(x = "Time (since exposure)", 
       y = "Ct value", 
       title = "Covariate-level prior predictive distribution") + 
  theme_minimal()

# saving the plot
ggsave("outputs/figures/pdfs/supplement/pop_prior_predictive.pdf", 
       p_prior_predictive, 
       width = 8,
       height = 5,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/pop_prior_predictive.png", 
       p_prior_predictive, 
       width = 8,
       height = 5,
       bg = "white")

