library(tidybayes)
library(ggridges)

source("scripts/setup/main.R")

n_obs <- obs[, uniqueN(data_id)]

# need to make plot of individual-level priors
dt_pop_priors_wide <- sample_pop_priors(
  50000,
  c_lod = 40,
  switch = FALSE, 
  data_format = "wide", 
  scale_type = "model")

# number of samples set to equal the number of individuals included in the study
dt_pop_priors_long <- sample_pop_priors(
  10000,
  c_lod = 40,
  switch = FALSE, 
  data_format = "long", 
  scale_type = "model")

# adding a sample ID, so the prior samples are easier to match with posterior
# samples
dt_pop_priors_long[, sample_id := 1:.N, by = parameter]

# moving to long format and samples the individual-level variation parameters
dt_ind_priors_long <- dt_pop_priors_long[, 
                   .(value = rnorm(n_obs, value, 0.2) + 
                       rtruncnorm(n_obs, 0, 0.2, a = 0, b = Inf)),
                   by = c("sample_id", "parameter")
                   ][, id := 1:.N, by = c("sample_id", "parameter")]

# moving back to wide format
dt_ind_priors_wide <- dcast(dt_ind_priors_long,
                            id + sample_id ~ parameter,
                            value.var = "value") |> 
  transform_to_natural()

# moving back to long format and transforming back to the natural scale for 
# the timing and Ct value parameter samples
dt_ind_priors_long_nat <- melt(dt_ind_priors_wide, 
                               id.vars = c("id", "sample_id"),
                               variable.name = "parameter")

# sampling from the infection time priors
dt_t_inf <- obs[, .(
  value = rnorm(10000, max(-onset_time, 0, na.rm = TRUE) + 5, 5)),
  by = "id"
  ][, parameter := factor("t_inf")
  ][, sample_id := 1:.N, by = c("id", "parameter")]

# combining data.tables of prior samples
dt_ind_priors_long_nat <- rbind(dt_ind_priors_long_nat, dt_t_inf)

# changing IDs and names to factors
dt_ind_priors_long_nat[, id := factor(id)]
dt_ind_priors_long_nat[, sample_id := factor(sample_id)]
dt_ind_priors_long_nat[, parameter := factor(parameter)]

# moving to wide format
dt_ind_priors_wide_nat <- dcast(dt_ind_priors_long_nat,
                                id + sample_id ~ parameter,
                                value.var = "value")

# transforming the timing parameters to time since infection
dt_ind_priors_wide_nat[, t_p := t_p - t_inf]
dt_ind_priors_wide_nat[, t_lod := t_lod - t_inf]
dt_ind_priors_wide_nat[, t_inf_plot := -t_inf]

# moving bakc to long format
dt_ind_priors_long_nat <- melt(dt_ind_priors_wide_nat, 
     id.vars = c("id", "sample_id"),
     variable.name = "parameter")

# saving VOC data (by ID) to merge with individual-level prior data
onset_data <- obs[, onset_time, by = c("id")] %>% unique()

voc_data <- obs[, VOC, by = c("id")] %>%
  unique()

# merging prior samples with VOC data from full dataset
dt_ind_priors_long_nat <- merge(dt_ind_priors_long_nat, 
                                voc_data[, id := factor(id)], 
                                by = "id")

# making sure the names of all variables match, to we can merge easily
dt_ind_priors_plot <- setnames(
  dt_ind_priors_long_nat, "parameter", "variable")

# changing the names of the parameters so that they match
dt_ind_priors_plot[, type := "prior"]
dt_ind_priors_plot[variable == "t_inf_plot", variable := "Time of infection"]
dt_ind_priors_plot[variable == "t_p", variable := "Time of peak"]
dt_ind_priors_plot[variable == "t_lod", variable := "Time LOD reached"]

#--- loading the main fit object to plot the individual-level posteriors
fit_main <- readRDS("outputs/fits/main.rds")

# extracting the individual-level posterior draws from the fit object
dt_ind_wide <- spread_draws(
  fit_main$draws(),
  inf_rel[id],
  c_0,
  t_inf[id],
  t_p[id],
  t_lod[id],
  c_p[id],
  sigma) %>% 
  data.table()

# adding c_lod, which we assume is the same as c_0 
dt_ind_wide[, c_lod := c_0]

# merging with onset data
dt_ind_wide <- merge(
  dt_ind_wide[, id := factor(id)], onset_data, by = "id")

# adjusting onset dates so they are relative to inferred infection times
dt_ind_wide_adj <- copy(dt_ind_wide)
dt_ind_wide_adj[, onset_time_abs := onset_time - quantile(t_inf, 0.5),
                by = id]

# make timing of peak relative to first positive test rather than relative to
# the inferred exposure time, which is what the posterior estimates are 
# relative to as they come out of the inference
dt_ind_wide_adj[, t_p := t_p - t_inf]
dt_ind_wide_adj[, t_lod := t_lod - t_inf]
dt_ind_wide_adj[, t_inf_plot := -t_inf]

# calculating the difference between the symptom onset (relative to first
# positive test) and the timing of the peak for each individual
dt_ind_wide_adj[, diff_t_p_onset := abs(onset_time - t_p), by = id]

# metling for plotting
dt_ind_long_adj <- melt(
  dt_ind_wide_adj,
  id.vars = c("id", ".chain", ".iteration", ".draw"))

# merging with VOC data 
dt_ind_long_adj <- merge(dt_ind_long_adj, voc_data, by = "id")

# timing posteriors plot
dt_t_plot <- dt_ind_long_adj[
  variable %in% c("t_inf_plot",
                  "t_p",
                  "t_lod",
                  "onset_time")
][, VOC := factor(
  VOC, 
  levels = c("Delta", 
             "Omicron (BA.1)", 
             "Omicron (BA.2)"),
  labels = c("Delta", "BA.1", "BA.2"))
][, VOC := fct_relevel(
  VOC, 
  c("Delta",
    "BA.1",
    "BA.2"))
][, variable := factor(
  variable,
  levels = c("t_inf_plot",
             "t_p",
             "t_lod",
             "onset_time"),
  labels = c("Time of infection",
             "Time of peak",
             "Time LOD reached",
             "Symptom onset"))]

# adding posterior type and removing redundant columns
dt_t_plot_final <- dt_t_plot[, type := "posterior"]
dt_t_plot_final[, .chain := NULL]
dt_t_plot_final[, .iteration := NULL]
dt_t_plot_final[, .draw := NULL]
dt_t_plot_final[, sample_id := 1:.N, by = c("id", "variable")]

# Making the names of the VOCs shorter for ease with plotting
dt_ind_priors_plot[VOC == "Omicron (BA.1)", VOC := "BA.1"]
dt_ind_priors_plot[VOC == "Omicron (BA.2)", VOC := "BA.2"]

# combining the prior samples with the posterior samples
dt_ind_prior_plot_all <- rbind(dt_t_plot_final, dt_ind_priors_plot)

# Just keeping the variables of interest for this figure
dt_t_prior_plot_all <- dt_ind_prior_plot_all[
  variable %in% c("Time of infection",
                  "Time of peak",
                  "Time LOD reached",
                  "Symptom onset")]

# isolating the IDs
# ids_ba1 <- dt_t_prior_plot_all[VOC == "BA.1", unique(id)]
# ids_plot <- sample(ids_ba1, 16)

# Plot for Delta individuals
p_delta <- ggplot() + 
  geom_density(data = dt_t_prior_plot_all[VOC == "Delta" & 
                                  variable != "Symptom onset"],
               aes(x = value, colour = variable, fill = type), alpha = 0.5) + 
  geom_vline(data = dt_t_prior_plot_all[VOC == "Delta" &
                                variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 25)) +
  labs(x = "Time relative to first positive test",
       y = "Density",
       title = "Delta timing posterior distributions") 

# saving Delta plot
ggsave("outputs/figures/pdfs/ind_prior_vs_posterior_delta.pdf", 
       p_delta,
       width = 10,
       height = 10)

# Plot for BA.1 individuals
p_ba1 <- ggplot() + 
  geom_density(data = dt_t_prior_plot_all[VOC == "BA.1" & 
                                            variable != "Symptom onset"],
               aes(x = value, colour = variable, fill = type), alpha = 0.5) + 
  geom_vline(data = dt_t_prior_plot_all[VOC == "BA.1" &
                                          variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 25)) +
  labs(x = "Time relative to first positive test",
       y = "Density",
       title = "Omicron (BA.1) timing posterior distributions") 

# saving BA.1 plot
ggsave("outputs/figures/pdfs/ind_prior_vs_posterior_ba1.pdf", 
       p_ba1,
       width = 10,
       height = 10)

# Plot for Delta individuals
p_ba2 <- ggplot() + 
  geom_density(data = dt_t_prior_plot_all[VOC == "BA.2" & 
                                          variable != "Symptom onset"],
               aes(x = value, colour = variable, fill = type), alpha = 0.5) + 
  geom_vline(data = dt_t_prior_plot_all[VOC == "BA.2" & 
                                        variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 25)) +
  labs(x = "Time relative to first positive test",
       y = "Density",
       title = "Omicron (BA.2) timing posterior distributions") 

# saving BA.2 plot
ggsave("outputs/figures/pdfs/ind_prior_vs_posterior_ba2.pdf", 
       p_ba2,
       width = 10,
       height = 10)


