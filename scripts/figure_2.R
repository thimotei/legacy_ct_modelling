library(lemon)
library(forcats)
library(cowplot)

source("scripts/setup.R")

# figure 2, panel A
dt_obs <- obs[no_pos_results >= 2]
dt_plot <- reshape_figure_2_panel_a(dt_obs)

p1_fig2 <- figure_2_panel_a(dt_plot) + 
  theme(legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin()) +
  labs(tag = "A") 

# figure 2, panel B
p2_fig2 <- figure_2_panel_b(obs[,
                                VOC := fct_relevel(VOC, "Delta")
], voc_arg = c("Delta",
               "Omicron (BA.1)",
               "Omicron (BA.2)")
) + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin()) +
  labs(tag = "B")

fig2 <- plot_grid(p1_fig2, p2_fig2)

ggsave("outputs/figures/pngs/figure_2.png",
       fig2,
       width = 12,
       height = 12,
       bg = "white")

ggsave("outputs/figures/pdfs/figure_2.pdf",
       fig2,
       width = 12,
       height = 12,
       bg = "white")
