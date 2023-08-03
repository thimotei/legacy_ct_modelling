# need to make plot of individual-level priors
dt_pop_priors_wide <- sample_pop_priors(50000,
                                        c_lod = 40,
                                        switch = FALSE, 
                                        data_format = "wide", 
                                        scale_type = "model")

# number of samples set to equal the number of individuals included in the study
dt_pop_priors_long <- sample_pop_priors(10000,
                                        c_lod = 40,
                                        switch = FALSE, 
                                        data_format = "long", 
                                        scale_type = "model")

dt_pop_priors_long[, sample_id := 1:.N, by = parameter]

dt_ind_priors_long <- dt_pop_priors_long[, 
                   .(value = rnorm(157, value, 0.2) + 
                       rtruncnorm(157, 0, 0.2, a = 0, b = Inf)),
                   by = c("sample_id", "parameter")
                   ][, id := 1:.N, by = c("sample_id", "parameter")]

dt_ind_priors_wide <- dcast(dt_ind_priors_long,
                            id + sample_id ~ parameter,
                            value.var = "value") |> 
  transform_to_natural()

dt_ind_priors_long_nat <- melt(dt_ind_priors_wide, 
                               id.vars = c("id", "sample_id"),
                               variable.name = "parameter")

dt_t_inf <- obs[, .(
  value = rnorm(10000, max(-onset_time, 0, na.rm = TRUE) + 5, 5)),
  by = "id"
  ][, parameter := factor("t_inf")
  ][, sample_id := 1:.N, by = c("id", "parameter")]

dt_ind_priors_long_nat <- rbind(dt_ind_priors_long_nat, dt_t_inf)

dt_ind_priors_long_nat[, id := factor(id)]
dt_ind_priors_long_nat[, sample_id := factor(sample_id)]
dt_ind_priors_long_nat[, parameter := factor(parameter)]

dt_ind_priors_wide_nat <- dcast(dt_ind_priors_long_nat,
                                id + sample_id ~ parameter,
                                value.var = "value")

dt_ind_priors_wide_nat[, t_p := t_p - t_inf]
dt_ind_priors_wide_nat[, t_lod := t_lod - t_inf]
dt_ind_priors_wide_nat[, t_inf_plot := -t_inf]

dt_ind_priors_long_nat <- melt(dt_ind_priors_wide_nat, 
     id.vars = c("id", "sample_id"),
     variable.name = "parameter")

dt_ind_priors_long_nat <- merge(dt_ind_priors_long_nat, 
                                voc_data[, id := factor(id)], 
                                by = "id")

dt_ind_priors_plot <- setnames(
  dt_ind_priors_long_nat, "parameter", "variable")

dt_ind_priors_plot[, type := "prior"]
dt_ind_priors_plot[variable == "t_inf_plot", variable := "Time of infection"]
dt_ind_priors_plot[variable == "t_p", variable := "Time of peak"]
dt_ind_priors_plot[variable == "t_lod", variable := "Time LOD reached"]

dt_t_plot_final <- dt_t_plot[, type := "posterior"]
dt_t_plot_final[, .chain := NULL]
dt_t_plot_final[, .iteration := NULL]
dt_t_plot_final[, .draw := NULL]
dt_t_plot_final[, sample_id := 1:.N, by = c("id", "variable")]

dt_ind_priors_plot[VOC == "Omicron (BA.1)", VOC := "BA.1"]
dt_ind_priors_plot[VOC == "Omicron (BA.2)", VOC := "BA.2"]

dt_ind_prior_plot_all <- rbind(dt_t_plot_final, dt_ind_priors_plot)

dt_t_prior_plot_all <- dt_ind_prior_plot_all[
  variable %in% c("Time of infection",
                  "Time of peak",
                  "Time LOD reached",
                  "Symptom onset")]

ggplot() + 
  geom_density(data = dt_t_prior_plot_all[id %in% 1:40 & VOC == "BA.1" & 
                                  variable != "Symptom onset"],
               aes(x = value, colour = variable, fill = type), alpha = 0.5) + 
  geom_vline(data = dt_t_prior_plot_all[id %in% 1:40 & VOC == "BA.1" &
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

dt_t_prior_plot_all[VOC == "Omicron (BA.1)"]


