figure_gene_targets <- function(dt, voc_arg){
  
  obs <- data.table::copy(dt)
  
  obs_plot <- copy(obs)[
    , no_vaccines := forcats::fct_relevel(no_vaccines, "2")
  ][
    , VOC := forcats::fct_relevel(VOC, "Delta")
  ]
  
  obs <- obs[VOC %in% voc_arg]
  
  obs[, combo.id := paste0(id,".",infection_id)]
  obs[, clean.id := .GRP, by = "combo.id"][, clean.id := as.factor(clean.id)]
  obs[, last_swab := max(t_since_first_pos), by = c("clean.id")]
  
  plot_data <- obs[t_since_first_pos >= -10]
  plot_data[, clean.id := forcats::fct_reorder(clean.id, last_swab)]
  
  p <- plot_data %>%
    ggplot(aes(x = t_since_first_pos,
               y = clean.id,
               col = ct_type)) +
    geom_line() +
    geom_point(aes(shape = result), alpha = 0.5, size = 2.5) +
    scale_shape_manual(values = c(15, 16), name = "") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(-10, 30, 5),
                       limits = c(NA, 30)) +
    scale_y_discrete(labels = paste0(plot_data$no_vaccines, ",",
                                     plot_data$total_infections, ",",
                                     plot_data$no_exposures)) + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    labs(x = "Days since first positive test") +
    facet_rep_wrap(~VOC, scales = "free")
  
  return(p)
}
