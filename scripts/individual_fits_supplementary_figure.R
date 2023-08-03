library(tidybayes)
library(data.table)
library(ggplot2)
library(forcats)
library(ggridges)

source("scripts/setup.R")

# load object with all fitted draws
fit <- readRDS("outputs/fits/fit_main.rds")

fit <- fit_without_voc

ct_model <- ct_model_no_voc

# subsetting onset and VOC data for merging later
onset_data <- obs[, onset_time, by = c("id")] %>% unique()
voc_data <- obs[, VOC, by = c("id")] %>% unique()

# extracting the individual-level posterior draws from the fit object
dt_ind_wide <- spread_draws(fit$draws(),
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
dt_ind_wide <- merge(dt_ind_wide, onset_data, by = "id")

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
dt_ind_long_adj <- melt(dt_ind_wide_adj,
     id.vars = c("id", ".chain", ".iteration", ".draw"))

# merging with VOC data 
dt_ind_long_adj <- merge(dt_ind_long_adj, voc_data, by = "id")

# timing posteriors plot
dt_t_plot <- dt_ind_long_adj[variable %in% c("t_inf_plot",
                                       "t_p",
                                       "t_lod",
                                       "onset_time")][,
  VOC := factor(VOC, 
                levels = c("Delta", "Omicron (BA.1)", "Omicron (BA.2)"),
                labels = c("Delta", "BA.1", "BA.2"))][,
  VOC := fct_relevel(VOC, c("Delta",
                            "BA.1",
                            "BA.2"))][,
  variable := factor(variable,
                     levels = c("t_inf_plot",
                                "t_p",
                                "t_lod",
                                "onset_time"),
                     labels = c("Time of infection",
                                "Time of peak",
                                "Time LOD reached",
                                "Symptom onset"))]

# plot of all timing posteriors
p_t_delta <- ggplot() + 
  geom_density(data = dt_t_plot[VOC == "Delta" & 
                                variable != "Symptom onset"],
               aes(x = value, fill = variable), alpha = 0.5) + 
  geom_vline(data = dt_t_plot[VOC == "Delta" & 
                              variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 30)) +
  labs(x = "Time relative to first positive test",
       y = "Density",
       title = "Delta timing posterior distributions") 

p_t_ba1 <- ggplot() + 
  geom_density(data = dt_t_plot[VOC == "BA.1" & 
                                  variable != "Symptom onset"],
               aes(x = value, fill = variable), alpha = 0.5) + 
  geom_vline(data = dt_t_plot[VOC == "BA.1" & 
                                variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 30)) +
  labs(x = "Time relative to first positive test", 
       y = "Density",
       title = "Omicron (BA.1) timing posterior distributions") 

p_t_ba2 <- ggplot() + 
  geom_density(data = dt_t_plot[VOC == "BA.2" & 
                                  variable != "Symptom onset"],
               aes(x = value, fill = variable), alpha = 0.5) + 
  geom_vline(data = dt_t_plot[VOC == "BA.2" & 
                                variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(id), scales = "free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-15, 30)) +
  labs(x = "Time relative to first positive test", 
       y = "Density",
       title = "Omicron (BA.2) timing posterior distributions") 

ggsave("outputs/figures/pdfs/delta_t_posteriors.pdf",
       p_t_delta,
       width = 10,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/ba1_t_posteriors.pdf",
       p_t_ba1,
       width = 10,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/ba2_t_posteriors.pdf",
       p_t_ba2,
       width = 10,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pngs/delta_t_posteriors.png",
       p_t_delta,
       width = 10,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pngs/ba1_t_posteriors.png",
       p_t_ba1,
       width = 10,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pngs/ba2_t_posteriors.png",
       p_t_ba2,
       width = 10,
       height = 10,
       bg = "white")

# Ct value posteriors plot
dt_ct_plot <- dt_ind_long_adj[variable %in% c("c_p")][, 
  VOC := factor(VOC, 
                levels = c("Delta", "Omicron (BA.1)", "Omicron (BA.2)"),
                labels = c("Delta", "BA.1", "BA.2"))][,
  VOC := fct_relevel(VOC, c("Delta",
                            "BA.1",
                            "BA.2"))][,
  variable := factor(variable,
                     levels = "c_p",
                     labels = "Ct value at peak")]

dt_ct_plot[, id := factor(id)]
dt_ct_plot[, id := fct_reorder(id, value, .desc = TRUE)]

p_ct <- dt_ct_plot %>% 
  ggplot() + 
  geom_density_ridges(aes(x = value, 
                          y = factor(id),
                          fill = VOC), 
                      alpha = 0.6) + 
  facet_wrap(vars(VOC), scales = "free", ncol = 8) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  # lims(x = c(10, 25)) +
  labs(x = "Ct value at peak",
       y = "Probability",
       title = "Ct value posterior distributions by VOC") + 
  scale_fill_brewer(palette = "Set1", aesthetics = "fill")

ggsave("outputs/figures/pdfs/ct_posteriors.pdf",
       p_ct,
       width = 8,
       height = 12,
       bg = "white")

ggsave("outputs/figures/pdfs/ct_posteriors.pdf",
       p_ct,
       width = 8,
       height = 12,
       bg = "white")


# simulating Ct trajectories from individual-level posteriors
dt_sims <- simulate_cts(params = dt_ind_wide,
                        time_range = seq(0, 30, 1),
                        obs_noise = FALSE)

# summarising simulated trajectories
dt_sims_sum <- summarise_ct_traj(dt_sims, pop_flag = FALSE)
dt_sims_sum_all <- merge(dt_sims_sum, voc_data, by = "id")

# calculating maximum infection time for adjustment
dt_t_inf <- dt_ind_wide[, .(t_inf_med = quantile(t_inf, 0.5)), by = id]

# merging with median infection time
obs_adj <- merge(obs, dt_t_inf, by = "id")

# adjusting raw data so its on the same scale as inferred trajectories
obs_adj[, t_first_test_since_inf := t_first_test + t_inf_med, by = "id"]
obs_adj[, onset_time_adj := onset_time + t_inf_med, by = "id"]

# relabelling factors for plot
obs_plot <- obs_adj[, ct_type := factor(ct_type,
                                    labels = c("ORF1ab",
                                               "N gene",
                                               "S gene"))][,
  VOC := factor(VOC, 
                levels = c("Delta", "Omicron (BA.1)", "Omicron (BA.2)"),
                labels = c("Delta", "BA.1", "BA.2"))][,
  VOC := fct_relevel(VOC, c("Delta", "BA.1", "BA.2"))]

dt_sims_sum_all[, 
  VOC := factor(VOC, 
                levels = c("Delta", "Omicron (BA.1)", "Omicron (BA.2)"),
                labels = c("Delta", "BA.1", "BA.2"))][,
  VOC := fct_relevel(VOC, c("Delta", "BA.1", "BA.2"))]

# adding posterior predictions
pp_summary <- summarise_pp(fit, obs_adj)

p_delta_fits <- 
  ggplot() + 
  geom_line(data = dt_sims_sum_all[VOC == "Delta"], 
            aes(x = t, y = median, group = id), 
            alpha = 0.5, linetype = "dashed") +
  geom_ribbon(data = dt_sims_sum_all[VOC == "Delta"],
              aes(x = t,
                  ymin = lo90,
                  ymax = hi90,
                  fill = VOC), 
              alpha = 0.5) +
  scale_fill_brewer(palette = "Set1", aesthetics = "fill") +
  geom_pointrange(data = pp_summary[VOC == "Delta"],
                  aes(x = t_first_test_since_inf,
                      y = median,
                      ymin = lo90,
                      ymax = hi90,
                      colour = ct_type)) +
  geom_vline(data = obs_plot[VOC == "Delta"],
             aes(xintercept = onset_time_adj), linetype = "dashed") +
  facet_wrap(vars(id)) + 
  scale_y_reverse() +
  coord_cartesian(clip = "off", ylim = c(40, 10)) + 
  scale_colour_brewer(palette = "Set2", aesthetics = "colour") +
  labs(x = "Time since exposure",
       y = "Ct value",
       fill = "VOC",
       colour = "Target gene",
       title = "Model fits for Delta-infected individuals") +
  theme_minimal() + 
  theme(legend.position = "bottom")

p_ba1_fits <- 
  ggplot() + 
  geom_line(data = dt_sims_sum_all[VOC == "BA.1"], 
            aes(x = t, y = median, group = id), 
            alpha = 0.5, linetype = "dashed") +
  geom_ribbon(data = dt_sims_sum_all[VOC == "BA.1"],
              aes(x = t,
                  ymin = lo90,
                  ymax = hi90,
                  fill = VOC), 
              alpha = 0.5) +
  scale_fill_brewer(palette = "Set1", aesthetics = "fill") +
  geom_pointrange(data = pp_summary[VOC == "BA.1"],
                  aes(x = t_first_test_since_inf,
                      y = median,
                      ymin = lo90,
                      ymax = hi90,
                      colour = ct_type)) +
  geom_vline(data = obs_plot[VOC == "BA.1"],
             aes(xintercept = onset_time_adj), linetype = "dashed") +
  facet_wrap(vars(id)) + 
  scale_y_reverse() +
  coord_cartesian(clip = "off", ylim = c(40, 10)) + 
  scale_colour_brewer(palette = "Set2", aesthetics = "colour") +
  labs(x = "Time since exposure",
       y = "Ct value",
       fill = "VOC",
       colour = "Target gene",
       title = "Model fits for BA.1-infected individuals") +
  theme_minimal() + 
  theme(legend.position = "bottom")


p_ba2_fits <- 
  ggplot() + 
  geom_line(data = dt_sims_sum_all[VOC == "BA.2"], 
            aes(x = t, y = median, group = id), 
            alpha = 0.5, linetype = "dashed") +
  geom_ribbon(data = dt_sims_sum_all[VOC == "BA.2"],
              aes(x = t,
                  ymin = lo90,
                  ymax = hi90,
                  fill = VOC), 
              alpha = 0.5) +
  scale_fill_brewer(palette = "Set1", aesthetics = "fill") +
  geom_pointrange(data = pp_summary[VOC == "BA.2"],
                  aes(x = t_first_test_since_inf,
                      y = median,
                      ymin = lo90,
                      ymax = hi90,
                      colour = ct_type)) +
  geom_vline(data = obs_plot[VOC == "BA.2"],
             aes(xintercept = onset_time_adj), linetype = "dashed") +
  facet_wrap(vars(id)) + 
  scale_y_reverse() +
  coord_cartesian(clip = "off", ylim = c(40, 10)) + 
  scale_colour_brewer(palette = "Set2", aesthetics = "colour") +
  labs(x = "Time since exposure",
       y = "Ct value",
       fill = "VOC",
       colour = "Target gene",
       title = "Model fits for BA.2-infected individuals") +
  theme_minimal() + 
  theme(legend.position = "bottom")

  
# saving PDF plot
ggsave("outputs/figures/pdfs/delta_fits.pdf",
       p_delta_fits,
       width = 8,
       height = 8)

# saving PNG plot
ggsave("outputs/figures/pngs/delta_fits.png",
       p_delta_fits,
       width = 8,
       height = 8)

# saving PDF plot
ggsave("outputs/figures/pdfs/ba1_fits.pdf",
       p_ba1_fits,
       width = 8,
       height = 10)

# saving PNG plot
ggsave("outputs/figures/pngs/ba1_fits.png",
       p_ba1_fits,
       width = 8,
       height = 10)

# saving PDF plot
ggsave("outputs/figures/pdfs/ba2_fits.pdf",
       p_ba2_fits,
       width = 8,
       height = 10)

# saving PNG plot
ggsave("outputs/figures/pngs/ba2_fits.png",
       p_ba2_fits,
       width = 8,
       height = 10)

