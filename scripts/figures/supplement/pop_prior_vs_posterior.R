library(ggh4x)

source("scripts/setup/main.R")

# load object with all fitted draws
fit_main <- readRDS("outputs/fits/fit_main.rds")

# extracting draws
draws <- extract_draws(fit_main)

tidybayes::spread_draws(fit_main, t_lod_mean)

# adjusting draws with the estimated covariate posteriors
adj_draws <- adjust_params(
  draws, 
  design = ct_model$design, 
  onsets_flag = FALSE) 

# updating predictor labels so they are more readable
adj_draws <- adj_draws %>% update_predictor_labels()

# adding the baseline predictor values and names to the draws data.table
adj_draws <- add_baseline_to_draws(
  adj_draws, "baseline", onsets_flag = FALSE) |> 
  transform_to_model()

# moving to long format and adding a column to distinguish between prior and
# posterior draws
adj_draws_long <- melt(
  adj_draws[, c("c_0", "c_p", "t_p", "t_lod", "predictor")], 
  measure.vars = c("c_0", "c_p",
                   "t_p", "t_lod"),
  variable.name = "parameter")[, type := "Posterior"][!is.na(predictor)]

# extracting the names of the predictors, for ease with mapping the prior
# samples over each predictor
predictors <- adj_draws_long[
  !is.na(predictor), predictor] |> 
  unique()

# sampling from the population-level priors
n_samp <- 100000

dt_pop_priors_long <- sample_pop_priors(
  n_samples = n_samp,
  c_lod = 40,
  data_format = "long",
  scale_type = "natural")

# removing parameters for model with switch, which we don't use for the main
# set of models
dt_pop_priors_long[, t_s := NULL]
dt_pop_priors_long[, c_s := NULL]
dt_pop_priors_long[, type := "Prior"]

# adding the set of predictors as a column for ease with plotting and comparing
# with posterior samples
n_params <- dt_pop_priors_long[, length(unique(parameter))]
dt_pop_priors_long[
  , predictor := rep(
    predictors, n_samp*n_params/length(predictors))
  ][order(predictor, parameter)][, type := "Prior"]

# merging prior samples with posteriors
dt_pop_prior_comparison <- rbind(
  dt_pop_priors_long, adj_draws_long)

# updating parameter names for plot
dt_pop_prior_comparison[
  parameter == "c_p", parameter := "Ct value at peak"]

dt_pop_prior_comparison[
  parameter == "t_p", parameter := "Timing of peak"]

dt_pop_prior_comparison[
  parameter == "t_lod", parameter := "Timing LOD reached"]

# plotting priors and posteriors in the same panel for comparison, stratified
# by predictor and parameter
p_prior_vs_post <- dt_pop_prior_comparison[
  parameter %in% c("Ct value at peak",
                   "Timing of peak",
                   "Timing LOD reached")
  ][, type := fct_relevel(type, "Prior")] %>% 
  ggplot() + 
  geom_density(aes(x = value, fill = type), alpha = 0.5) + 
  lims(x = c(0, 40)) +
  labs(x = "Value (Ct value or Time (days))", y = "Density", fill = "Type") +
  facet_grid(parameter~predictor, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "bottom") 

# saving plots
ggsave("outputs/figures/pdfs/supplement/pop_prior_vs_posterior.pdf",
       p_prior_vs_post,
       width = 9,
       height = 6)

ggsave("outputs/figures/pngs/supplement/pop_prior_vs_posterior.png",
       p_prior_vs_post,
       width = 9,
       height = 6)
