library(tidybayes)
library(data.table)
library(ggplot2)
library(forcats)

# load object with all fitted draws
fit <- readRDS("outputs/fits/fit_full.rds")

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

# make timing of peak relative to first positive test
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
  VOC := fct_relevel(VOC, c("Delta",
                            "Omicron (BA.1)",
                            "Omicron (BA.2)"))][,
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
p_t_all <- ggplot() + 
  geom_density(data = dt_t_plot[variable != "Symptom onset"],
               aes(x = value, fill = variable), alpha = 0.5) + 
  geom_vline(data = dt_t_plot[variable == "Symptom onset"],
             aes(xintercept = value), linetype = "dashed") +
  facet_wrap(vars(VOC, id), scales = "free", ncol = 7) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(-10, 20)) +
  labs(x = "Time relative to first positive test", y = "Density")

ggsave("outputs/figures/pdfs/individual_t_posteriors.pdf",
       p_t_all,
       width = 10,
       height = 30,
       bg = "white")

ggsave("outputs/figures/pngs/individual_t_posteriors.png",
       p_t_all,
       width = 10,
       height = 30,
       bg = "white")

# Ct value posteriors plot
dt_ct_plot <- dt_ind_long_adj[variable %in% c("c_p")][,
  VOC := fct_relevel(VOC, c("Delta",
                            "Omicron (BA.1)",
                            "Omicron (BA.2)"))][,
  variable := factor(variable,
                     levels = "c_p",
                     labels = "Ct value at peak")]

p_ct_all <- dt_ct_plot %>% 
  ggplot(aes(x = value, fill = variable)) + 
  geom_density() + 
  facet_wrap(vars(VOC, id), scales = "free", ncol = 7) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  lims(x = c(10, 25), y = c(0, 0.5)) +
  labs(x = "Ct value at peak", y = "Probability density")

ggsave("outputs/figures/pdfs/individual_ct_posteriors.pdf",
       p_ct_all,
       width = 10,
       height = 30,
       bg = "white")

ggsave("outputs/figures/pngs/individual_ct_posteriors.png",
       p_ct_all,
       width = 10,
       height = 30,
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
                                               "S gene"))]

# plotting all individual-level inferred trajectories with raw data on top
p_all_fits <- dt_sims_sum_all[, VOC := fct_relevel(VOC,
                                     c("Delta",
                                       "Omicron (BA.1)",
                                       "Omicron (BA.2)"))] %>% 
 ggplot() + 
  geom_line(aes(x = t, y = median, group = id), 
            alpha = 0.5, linetype = "dashed") +
  geom_ribbon(aes(x = t, ymin = lo90, ymax = hi90, fill = VOC), 
              alpha = 0.5) +
  scale_fill_brewer(palette = "Set1", aesthetics = "fill") +
  geom_point(data = obs_plot,
             inherit.aes = FALSE,
             aes(x = t_first_test_since_inf, y = ct_value, colour = ct_type)) +
  geom_vline(data = obs_plot,
             aes(xintercept = onset_time_adj), linetype = "dashed") +
  facet_wrap(vars(VOC, id), 
             ncol = 7) + 
  scale_y_reverse() +
  coord_cartesian(clip = "off", ylim = c(40, 10)) + 
  scale_colour_brewer(palette = "Set2", aesthetics = "colour") +
  labs(x = "Time since exposure",
       y = "Ct value",
       fill = "VOC",
       colour = "Target gene") +
  theme_minimal() + 
  theme(legend.position = "bottom")
  
# saving PDF plot
ggsave("outputs/figures/pdfs/all_individuals_fits.pdf",
       p_all_fits,
       width = 10,
       height = 30)

# saving PNG plot
ggsave("outputs/figures/pngs/all_individuals_fits.png",
       p_all_fits,
       width = 10,
       height = 30)

