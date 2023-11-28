library(lemon)
library(forcats)
library(cowplot)
library(ggsci)

source("scripts/setup/main.R")

# changing factor levels
obs[, VOC := fct_relevel(VOC, "Delta")]
obs[, `Symptoms` := fct_relevel(symptoms, "symptomatic")]
obs[, `Number of exposures` := fct_relevel(no_exposures, "3")]
obs[, `Age group` := fct_relevel(age_group, "20-34")]

# figure 2, panel A
dt_obs <- obs[no_pos_results >= 2]
dt_plot <- reshape_figure_2_panel_a(dt_obs)

p_fig_2_a <- figure_2_panel_a(dt_plot) + 
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin()) +
  labs(tag = "A")

p_fig_2_b <- obs[, VOC := fct_relevel(VOC, "Delta")] |> 
  ggplot() +
  geom_point(aes(x = t_since_first_pos, y = ct_value, colour = VOC), alpha = 0.35) +
  geom_line(aes(x = t_since_first_pos, y = ct_value, colour = VOC, group = id), alpha = 0.045) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + 
  geom_hline(aes(yintercept = 37), linetype = "dashed") + 
  # geom_smooth(aes(x = t_since_first_pos, y = ct_value, colour = VOC), alpha = 0.1) +
  facet_wrap(~VOC) + 
  theme_minimal() + 
  scale_y_reverse() + 
  scale_colour_nejm() + 
  theme(legend.position = "none") + 
  labs(x = "Time (since first positive test)", 
       y = "Ct value", 
       title = "Observed Ct values by VOC")

p_fig_2_a_b <- plot_grid(p_fig_2_a, p_fig_2_b,
                         ncol = 1,
                         rel_heights = c(1, 1))

dt_obs_long <- melt(obs[, .(`id`, `ct_value`, `result`, `VOC`, 
                            `symptoms`, `no_exposures`, `age_group`)],
                    measure.vars = c("VOC", "symptoms", "no_exposures", "age_group"),
                    variable.name = "Regressor category",
                    value.name = "Regressor")

dt_obs_plot <- dt_obs_long[result == "Positive"]

p_fig_2_c_1 <- figure_2_ct_box_plot(dt_obs_plot, "VOC", "Delta") +
  labs(x = "") +
  scale_colour_nejm() +
  scale_fill_nejm() + 
  labs(title = "VOC")

p_fig_2_c_2 <- figure_2_ct_box_plot(dt_obs_plot, "symptoms", "symptomatic") +
  labs(x = "") +
  scale_colour_npg() +
  scale_fill_npg() +
  labs(title = "Symptom status")

p_fig_2_c_3 <- figure_2_ct_box_plot(dt_obs_plot, "age_group", "20-34") +
  labs(x = "") +
  scale_colour_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Age group")

p_fig_2_c_4 <- figure_2_ct_box_plot(dt_obs_plot, "no_exposures", "3") +
  labs(x = "") +
  scale_colour_aaas() +
  scale_fill_aaas() + 
  labs(title = "No. exposures")

p_fig_2_c <- plot_grid(p_fig_2_c_1, p_fig_2_c_2,
                       p_fig_2_c_3, p_fig_2_c_4)


p_fig_2_b_c <- plot_grid(p_fig_2_b, p_fig_2_c, ncol = 1,
                         rel_heights = c(0.66, 1))

p_fig_2 <- plot_grid(p_fig_2_a, p_fig_2_b_c, 
                     ncol = 2)

ggsave("outputs/figures/pngs/figure_2.png",
       p_fig_2,
       width = 14,
       height = 10,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_2.pdf",
       p_fig_2,
       width = 14,
       height = 10,
       bg = "white")

