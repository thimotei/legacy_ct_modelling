library(lemon)
library(forcats)
library(cowplot)

source("scripts/setup/main.R")

# figure 2, panel B
p_fig_ct_values <- figure_2_panel_b(
  obs[, VOC := fct_relevel(VOC, "Delta")],
  voc_arg = c("Delta",
               "Omicron (BA.1)",
               "Omicron (BA.2)")) + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin())

ggsave("outputs/figures/pngs/supplement/ind_testing_ct_values.png",
       p_fig_ct_values,
       # width = ,
       # height = 6,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/ind_testing_ct_values.pdf",
       p_fig_ct_values,
       width = 12,
       height = 8,
       bg = "white")
