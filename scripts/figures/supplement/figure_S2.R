# Script for making Figure S2
# Numbers of tests and gene targets detected for each individual

# Loading data
source("scripts/setup/main.R")

# Subsetting the data for the plot
dt_plot <- obs[, .(
  id, infection_id, t_since_first_pos, VOC, 
  no_vaccines, total_infections, no_exposures, ct_value, 
  result)]

# figure 2, panel B
p_fig_ct_values <- figure_2_panel_b(
  dt_plot[, VOC := fct_relevel(VOC, "Delta")],
  voc_arg = c(
    "Delta", "Omicron (BA.1)", "Omicron (BA.2)")) + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin())

# Saving the plotted data
saveRDS(dt_plot, "outputs/plot_data/supplement/figure_S2.rds")

ggsave("outputs/figures/pngs/supplement/figure_S2.png",
       p_fig_ct_values,
       # width = ,
       # height = 6,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/figure_S2.pdf",
       p_fig_ct_values,
       width = 12,
       height = 8,
       bg = "white")
