library(ggh4x)

# Loading data, model structure, etc.
source("scripts/setup/main.R")

# Load object with all fitted draws
fit_main <- readRDS("outputs/fits/main.rds")

# Extracting draws
draws <- extract_draws(fit_main)

# Adjusting draws with the estimated covariate posteriors
adj_draws <- adjust_params(
  draws, 
  design = ct_model$design, 
  onsets_flag = FALSE) 

# Updating predictor labels so they are more readable
adj_draws <- adj_draws %>% update_predictor_labels()

# Adding the baseline predictor values and names to the draws data.table
adj_draws <- add_baseline_to_draws(
  adj_draws, "baseline", onsets_flag = FALSE) |> 
  transform_to_model()

# Moving to long format and adding a column to distinguish between prior and
# posterior draws
adj_draws_long <- melt(
  adj_draws[, c("c_0", "c_p", "t_p", "t_lod", "predictor")], 
  measure.vars = c("c_0", "c_p",
                   "t_p", "t_lod"),
  variable.name = "parameter")[, type := "Posterior"][!is.na(predictor)]

# Extracting the names of the predictors, for ease with mapping the prior
# samples over each predictor
predictors <- adj_draws_long[
  !is.na(predictor), predictor] |> 
  unique()

# Sampling from the population-level priors
n_samp <- 100000
dt_pop_priors_long <- sample_pop_priors(
  n_samples = n_samp,
  c_lod = 40,
  data_format = "long",
  scale_type = "natural")

# Removing parameters for model with switch, which we don't use for the main
# set of models
dt_pop_priors_long[, t_s := NULL]
dt_pop_priors_long[, c_s := NULL]
dt_pop_priors_long[, type := "Prior"]

# Adding the set of predictors as a column for ease with plotting and comparing
# with posterior samples
n_params <- dt_pop_priors_long[, length(unique(parameter))]
dt_pop_priors_long[
  , predictor := rep(
    predictors, n_samp*n_params/length(predictors))
  ][order(predictor, parameter)][, type := "Prior"]

# Merging prior samples with posteriors
dt_pop_prior_comparison <- rbind(
  dt_pop_priors_long, adj_draws_long)

# updating parameter names for plot
dt_pop_prior_comparison[
  parameter == "c_p", parameter := "Ct value at peak"]

dt_pop_prior_comparison[
  parameter == "t_p", parameter := "Timing of peak"]

dt_pop_prior_comparison[
  parameter == "t_lod", parameter := "Timing LOD reached"]

# Plotting priors and posteriors in the same panel for comparison, stratified
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

# Saving data required to remake plot
saveRDS(dt_pop_prior_comparison, "outputs/plot_data/supplement/figure_S24.rds")

# Saving plots
ggsave("outputs/figures/pdfs/supplement/figure_S24.pdf",
       p_prior_vs_post,
       width = 9,
       height = 6,
       bg = "white")

ggsave("outputs/figures/pngs/supplement/figure_S24.png",
       p_prior_vs_post,
       width = 9,
       height = 6,
       bg = "white")
