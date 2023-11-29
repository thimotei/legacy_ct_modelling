# Loading data
source("scripts/setup/main.R")

# Plotting the number of observations available per patient
p1_timing <- obs[, .N, by = c("data_id", "VOC", "result")][
  , VOC := fct_relevel(VOC, "Delta")] |> 
  ggplot() +
  geom_bar(aes(x = N, fill = VOC)) +
  theme_minimal() + 
  labs(x = "Number of PCR tests per id", 
       y = "Count",
       tag = "A") + 
  facet_wrap(~result, scales = "free") +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Plotting the observed delay between first viral load data and symptom onset
p2_timing <- obs[, .(onset_time), by = c("data_id", "VOC")][
  , VOC := fct_relevel(VOC, "Delta")] |>
  unique() |> 
  ggplot() +
  geom_bar(aes(x = onset_time, fill = VOC)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + 
  theme_minimal() + 
  labs(title = "Time of symptom onset", 
       x = "Time (days since first +ve)",
       y = "Count",
       tag = "B") +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Plotting the observed peak viral load
p3_timing <- obs[, .SD[which.min(ct_value)], by = id
][, VOC := fct_relevel(VOC, "Delta")] |>
  ggplot(aes(x = VOC, y = ct_value, colour = VOC, fill = VOC)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.35) +
  geom_boxplot(aes(fill = VOC), 
               alpha = 0.2, outlier.shape = NA) + 
  scale_y_reverse() + 
  theme_minimal() + 
  labs(y = "Ct value",
       title = "Lowest observed Ct values",
       tag = "C") + 
  theme(legend.position = "none") +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Plotting the observed time to clearance
p4_timing <- obs[result == "Positive", .SD[which.max(t)], by = id
][, VOC := fct_relevel(VOC, "Delta")] |>
  ggplot(aes(x = VOC, y = t, colour = VOC, fill = VOC)) +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.35) + 
  geom_boxplot(aes(fill = VOC),
               alpha = 0.2, outlier.shape = NA) + 
  theme_minimal() + 
  labs(y = "Time (days since first +ve test)", 
       title = "Time of last +ve PCR test",
       tag = "D") + 
  theme(legend.position = "none") +
  scale_colour_nejm() + 
  scale_fill_nejm()

# Putting together the various panels for the final plot
p_top_row <- p1_timing + p2_timing + plot_layout(guides = "collect")
p_bottom_row <- p3_timing + p4_timing
p_timing <- plot_grid(p_top_row, p_bottom_row, nrow = 2)

# Saving final plot
ggsave("outputs/figures/pngs/supplement/figure_S5.png", 
       p_timing,
       width = 10,
       height = 7,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/figure_S5.pdf", 
       p_timing,
       width = 10,
       height = 7,
       bg = "white")
