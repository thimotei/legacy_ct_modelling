library(lemon)
library(forcats)
library(cowplot)

source("scripts/setup.R")

# changing factor levels
obs[, VOC := fct_relevel(VOC, "Delta")]
obs[, `Symptoms` := fct_relevel(symptoms, "symptomatic")]
obs[, `Number of exposures` := fct_relevel(no_exposures, "3")]
obs[, `Age group` := fct_relevel(age_group, "20-34")]

# figure 2, panel A
dt_obs <- obs[no_pos_results >= 2]
dt_plot <- reshape_figure_2_panel_a(dt_obs)

p_fig_2_a <- figure_2_panel_a(dt_plot, flip_arg = TRUE) + 
  theme(legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin()) +
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

p_fig_2_c_1 <- obs[, .SD[which.min(ct_value)], by = "id"] |> 
  unique() |> 
  ggplot(aes(x = VOC, y = ct_value, colour = VOC)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.35) +
  geom_boxplot(aes(fill = VOC), 
               alpha = 0.2, outlier.shape = NA) + 
  scale_y_reverse() + 
  theme_minimal() + 
  labs(y = "Ct value",
       title = "Lowest observed Ct values") + 
  theme(legend.position = "none") +
  scale_colour_nejm() + 
  scale_fill_nejm()

p_fig_2_c_2 <- obs[no_pos_results >= 2, .SD[which.max(t)], by = "id"] |> 
  ggplot(aes(x = VOC, y = t, colour = VOC)) +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.35) + 
  geom_boxplot(aes(fill = VOC),
               alpha = 0.2, outlier.shape = NA) + 
  theme_minimal() + 
  labs(y = "Time (days since first +ve test)", 
       title = "Time of last +ve PCR test") + 
  theme(legend.position = "none") +
  scale_colour_nejm() + 
  scale_fill_nejm()

p_fig_2_c_3 <- obs[onset_time > -15][, onset_time, by = c("id", "VOC")] |>
  unique() |> 
  ggplot() + 
  geom_bar(aes(x = onset_time, fill = VOC)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + 
  theme_minimal() + 
  labs(x = "Time (days since first +ve)", y = "",
       title = "Time of symptom onset") + 
  scale_colour_nejm() + 
  scale_fill_nejm() + 
  theme(legend.position = "none")

p_fig_2_c <- plot_grid(p_fig_2_c_1, p_fig_2_c_2, 
                       p_fig_2_c_3, nrow = 1)

p_fig_2_b_c <- plot_grid(p_fig_2_b, p_fig_2_c,
                         ncol = 1,
                         rel_heights = c(1, 1))

p_fig_2 <- plot_grid(p_fig_2_b_c, p_fig_2_a,
                     ncol = 1, rel_widths = c(1, 1))


ggsave("outputs/figures/pngs/figure_2_new_draft_2.png",
       p_fig_2,
       width = 12,
       height = 12,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_2_new_draft_2.pdf",
       p_fig_2,
       width = 12,
       height = 12,
       bg = "white")

