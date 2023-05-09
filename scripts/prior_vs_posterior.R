adj_draws <- adjust_params(draws, 
                           design = ct_model$design, 
                           onsets_flag = FALSE) 

adj_draws <- adj_draws %>% update_predictor_labels()

adj_draws <- add_baseline_to_draws(
  adj_draws, "baseline", onsets_flag = FALSE) %>% 
  transform_to_natural()

adj_draws_long <- melt(
  adj_draws[, c("c_0", "c_p", "t_p", "t_lod", "predictor")], 
  measure.vars = c("c_0", "c_p",
                   "t_p", "t_lod"),
  variable.name = "parameter")[, type := "Posterior"][!is.na(predictor)]


predictors <- adj_draws_long[!is.na(predictor), predictor] %>% unique()

dt_pop_prior_samples_long <- map_df(
  predictors, 
  sample_from_pop_priors, 
  samples = 100000,
  c_lod = 40)

dt_pop_prior_comparison <- rbind(dt_pop_prior_samples_long, adj_draws_long)  

p_prior_vs_post <- dt_pop_prior_comparison[parameter %in% c("c_p", "t_p", "t_lod")] %>% 
  ggplot() + 
  geom_density(aes(x = value, fill = type), alpha = 0.5) + 
  lims(x = c(0, 50)) +
  labs(x = "Value (Ct value or Time (days))", y = "Density", fill = "Type") +
  facet_grid(parameter~predictor, scales = "free") + 
  theme(legend.position = "bottom") + 
  theme_minimal() 

ggsave("outputs/figures/pdfs/pop_prior_vs_posterior.pdf",
       p_prior_vs_post,
       width = 14,
       height = 9)

ggsave("outputs/figures/pngs/pop_prior_vs_posterior.png",
       p_prior_vs_post,
       width = 14,
       height = 9)
