#--- custom plot theme
custom_plot_theme = list(
  theme_cowplot(11),
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.title = element_blank()),
  geom_hline(aes(yintercept = -Inf)),
  geom_vline(aes(xintercept = -Inf))
)