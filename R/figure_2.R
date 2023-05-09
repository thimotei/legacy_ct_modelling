figure_2_panel_a <- function(dt) {
  
  dt_plot <- data.table::copy(dt)
  
  p_out <- dt_plot %>% 
    ggplot() + 
    geom_linerange(aes(xmin = start_date, 
                       xmax = end_date,
                       y = factor(new_id),
                       colour = dose,
                       group = interaction(new_id, infection_id)),
                   size = 1.35,
                   alpha = 0.75) + 
    labs(colour = "Dose") + 
    ggsci::scale_colour_nejm() +  
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    new_scale_colour() +
    geom_point(aes(x = date,
                   y = new_id,
                   group = new_id,
                   colour = `VOC/event`,
                   shape = `VOC/event`),
               size = 2) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    # scale_x_continuous(breaks = c(ymd("2021-01-01"),
    #                               ymd("2021-07-01"),
    #                               ymd("2022-01-01")),
    #                    labels = c("January 2021",
    #                               "July 2021",
    #                               "Janurary 2022")) + 
    labs(x = "Date", y = "Participant") + 
    scale_colour_brewer(palette = "Set1") 
  
  return(p_out)
  
}

figure_2_panel_b <- function(dt, voc_arg){

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
               col = ct_value)) +
    geom_line(alpha = 0.5) +
    geom_point(aes(shape = result), size = 2.5) +
    scale_shape_manual(values = c(15, 16), name = "") +
    scale_colour_gradientn(
      colours = c("#377EB8", "#C21E56", "#FF0000"),
      breaks = c(10, 20, 30, 40),
      limits = c(40, 10),
      trans = c("reverse"),
      name = "Ct value") +
    # cowplot::theme_minimal_grid() +
    theme_minimal() +
    scale_x_continuous(breaks = seq(-10, 30, 5),
                       limits = c(NA, 30)) +
    scale_y_discrete(labels = paste0(plot_data$no_vaccines, ",",
                                     plot_data$total_infections, ",",
                                     plot_data$no_exposures)) + 
    # geom_text(aes(y = forcats::fct_reorder(clean.id, last_swab),
    #               x = -14,
    #               label = infection_id), inherit.aes = FALSE, size = 2.5) +
    # geom_text(aes(y = forcats::fct_reorder(clean.id, last_swab),
    #               x = 31,
    #               label = paste0(total_infections, ", ")), inherit.aes = FALSE, size = 2.5) +
    # geom_text(aes(y = forcats::fct_reorder(clean.id, last_swab),
    #               x = 32,
    #               label = no_vaccines), inherit.aes = FALSE, size = 2.5) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    labs(x = "Days since first positive test") +
    facet_rep_wrap(~VOC, scales = "free")

  return(p)
}
