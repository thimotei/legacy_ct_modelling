#--- custom plot theme
custom_plot_theme <- function(flip = FALSE) {
  
  custom_plot_theme = list(
    theme_cowplot(11),
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.title = element_blank()),
    geom_vline(aes(xintercept = -Inf)))

    if(flip == FALSE) { 
      append(custom_plot_theme, geom_hline(aes(yintercept = -Inf)))
    }
    
    return(custom_plot_theme)
}
