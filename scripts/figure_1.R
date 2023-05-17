# script for figure 1
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)

source("scripts/schematic_plot_a.R")
source("scripts/schematic_plot_b.R")

p1_fig1a <- schematic_plot_a() + 
  theme(text = element_text(family = "Fira Sans")) + 
  labs(tag = "")

p1_fig1b <- schematic_plot_b() + 
  theme(text = element_text(family = "Fira Sans"))

p1_figA <- p1_fig1a + p1_fig1b

ggsave("outputs/figures/fig1A.png",
       p1_figA,
       width = 8,
       height = 4,
       bg = "white")
