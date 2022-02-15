#--- custom plot theme
custom_plot_theme <- function(flip = FALSE, legend_arg = FALSE) {
  
  custom_plot_theme = list(
    theme_bw(11),
    theme(strip.placement = "outside"),
    geom_vline(aes(xintercept = -Inf)),
    guides(color = guide_legend(override.aes = list(fill = NA))))
  
  if(flip == FALSE) { 
    append(custom_plot_theme, geom_hline(aes(yintercept = -Inf)))
  }
  
  if(legend_arg == FALSE) {
    append(custom_plot_theme,
           theme(legend.title = element_blank()))
  }
  
  return(custom_plot_theme)
}
