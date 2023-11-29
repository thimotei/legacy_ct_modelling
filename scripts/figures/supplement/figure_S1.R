# Script for making Figure S1
# Numbers of tests and gene targets detected for each individual

# Loading data
source("scripts/setup/main.R")

# Subsetting the data for the plot
dt_plot <- obs[, .(
  id, infection_id, t_since_first_pos, VOC, 
  no_vaccines, total_infections, no_exposures, ct_value,
  ct_type, result)]

# Renaming the gene targets for plotting
dt_plot[ct_type == "ct_value", ct_type := "ORF1ab"]
dt_plot[ct_type == "ct_n_gene", ct_type := "S gene"]
dt_plot[ct_type == "ct_s_gene", ct_type := "N gene"]

# Plotting the data
p2_gene_targets <- figure_gene_targets(
  dt_plot[, VOC := fct_relevel(VOC, "Delta")],
  voc_arg = c("Delta",
               "Omicron (BA.1)",
               "Omicron (BA.2)")) + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin()) +
  labs(colour = "Gene target") + 
  scale_colour_brewer(palette = "Set1")

# Saving the plotted data
saveRDS(dt_plot, "outputs/plot_data/supplement/figure_S1.rds")

# Saving the plot
ggsave("outputs/figures/pngs/supplement/figure_S1.png",
       p2_gene_targets,
       width = 8,
       height = 8,
       bg = "white")

ggsave("outputs/figures/pdfs/supplement/figure_S1.pdf",
       p2_gene_targets,
       width = 8,
       height = 8,
       bg = "white")
