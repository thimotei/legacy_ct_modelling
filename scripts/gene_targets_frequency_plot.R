library(lemon)
library(forcats)
library(cowplot)

# figure 2, panel A
dt_plot <- obs[no_pos_results >= 2]

dt_plot[ct_type == "ct_value", ct_type := "ORF1ab"]
dt_plot[ct_type == "ct_n_gene", ct_type := "S gene"]
dt_plot[ct_type == "ct_s_gene", ct_type := "N gene"]

# figure 2, panel B
p2_figS1 <- figure_gene_targets(dt_plot[,
  VOC := fct_relevel(VOC, "Delta")
], voc_arg = c("Delta",
               "Omicron (BA.1)",
               "Omicron (BA.2)")
) + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin()) +
  labs(colour = "Gene target") + 
  scale_colour_brewer(palette = "Set1")

ggsave("outputs/figures/pdfs/figure_gene_targets_frequency.pdf",
       p2_figS1,
       width = 8,
       height = 8,
       bg = "white")

ggsave("outputs/figures/pngs/figure_gene_targets_frequency.png",
       p2_figS1,
       width = 8,
       height = 8,
       bg = "white")

