schematic_plot_a <- function() {
  
  # fake ct data and fits
  ct = data.table(time = c(0, 4.3, 18.5),
                  ct_value = c(40, 17.5, 40))
  
  ct.ci.lwr = splinefunH(c(1, 4.6, 17),
                         c(40, 24, 40),
                         c(0 ,0, 0))
  
  ct.ci.lwr = ct.ci.lwr(2:(17/0.01)*0.01)
  
  ct.ci.lwr = data.table(time = 2:(17/0.01)*0.01,
                         ct_low_cri = ct.ci.lwr)
  
  ct.ci.upp = splinefunH(c(-2, 4.5, 24),
                         c(40, 15, 40),
                         c(0, 0, 0))
  
  ct.ci.upp = ct.ci.upp(0:(20/0.01)*0.01)
  
  ct.ci.upp = data.table(time = 0:(20/0.01)*0.01,
                         ct_high_cri = ct.ci.upp)
  
  ct.ci.1 = merge(ct.ci.lwr, 
                  ct.ci.upp, all.y = TRUE)[is.na(ct_low_cri) == TRUE,
                                           ct_low_cri := 40][,
                                                             VOC := "Omicron"]
  
  ct.ci.2 = merge(ct.ci.lwr, 
                  ct.ci.upp, all.y = TRUE)
  
  ct.ci.2 = ct.ci.2[is.na(ct_low_cri) == TRUE,
                    ct_low_cri := 40][, VOC := "Delta"][,
  time := time - 1][,
  ct_low_cri := ct_low_cri + 0.5][, 
  ct_high_cri := ct_high_cri + 0.5]
  
  ct.ci = rbind(ct.ci.1, ct.ci.2)
  
  ct.obs = data.table(id = c(rep(1, 9), rep(2, 8)),
                      time = c(-6, -3, 0, 4, 6, 11, 15, 18, 21, 
                               -7, -4, 0, 3, 5, 9, 13, 16),
                      ct_value = c(40, 40, 36.5, 23, 20, 30, 35, 40, 40,
                                   40, 40, 35, 26.5, 18, 28, 34, 40),
                      type = rep("ct_value", 17),
                      # onset_time = c(4, 4, 4, 4, 4, 4, 4, 4, 4,
                      #           5, 5, 5, 5, 5, 5, 5, 5),
                      VOC = c(rep("Omicron", 9),
                              rep("Delta", 8)))
  
  ct_fit = data.table(id = c(1, 1, 1, 2, 2, 2),
                      time = c(0, 4.3, 18.5, 0, 5, 17),
                      ct_value = c(40, 17.5, 40, 40, 19, 40),
                      VOC = c(rep("Omicron", 3),
                              rep("Delta", 3)))
  # ct.obs = rbind(ct.obs, 
  #       data.table(id = 1, time = 4, ct_value = 40, type = "onset_time", VOC = "Delta"),
  #       data.table(id = 2, time = 5, ct_value = 40, type = "onset_time", VOC = "Omicron"))
  
  symptom_onset =  rbind(
    data.table(id = 1, time = 4, ct_value = 40, type = "onset_time", VOC = "Delta"),
    data.table(id = 2, time = 5.5, ct_value = 40, type = "onset_time", VOC = "Omicron")
    )
  
  # symptoms onset
  # symp.obs = data.table(time = 5, onset = 36)
  
  # colours
  # custom_pallete = list() 
  # CONFIG$cols = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF', '#F39B7FFF', 
  #                 '#8491B4FF', '#91D1C2FF', '#DC0000FF', '#7E6148FF', '#B09C85FF') # Nature
  
  # Replicating nature colour palette
  custom_pallete = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF', '#F39B7FFF', 
                     '#8491B4FF', '#91D1C2FF', '#DC0000FF', '#7E6148FF', '#B09C85FF')
  
  # changing colour to represent opacity for cartoon credible interval
  lightup = function(c, alpha)
  {
    z = col2rgb(c)/255
    return(rgb(z[1], z[2], z[3], alpha))  # 0.125 for col3
  }
  
  custom_pallete_ci_1= c()
  for(i in seq(custom_pallete)) {
    custom_pallete_ci_1[i] = lightup(custom_pallete[i], alpha = 0.25)
  }
  
  custom_pallete_ci_2 = c() 
  for(i in seq(custom_pallete)) {
    custom_pallete_ci_2[i] = lightup(custom_pallete[i], alpha = 0.1)
  }
  
  # sample plot
  figure_1_aa = ggplot(data = ct.obs) + 
    geom_point(aes(x = time, y = ct_value, colour = VOC), shape = 1) +
    annotate("segment", x = 4, y = 40, xend = 4, yend = 42.5,
             size = 0.5, linejoin = "mitre",
             arrow = arrow(type = "closed", length = unit(0.01, "npc"))) +
    annotate("segment", x = 6, y = 40, xend = 6, yend = 42.5,
             size = 0.5, linejoin = "mitre",
             arrow = arrow(type = "closed", length = unit(0.01, "npc"))) +
    # annotate("text", x = 4.5, y = 44, label = "symptom onset times",
    #          color = "black", angle = 0, hjust = 0.45, size = 2) +
    annotate("text", x = 4, y = 43.5, label = "Covariate 1 \n onset",
             color = "black", angle = 0, hjust = 0.45, size = 2) +
    annotate("text", x = 6, y = 43.5, label = "Covariate 2 \n onset",
             color = "black", angle = 0, hjust = 0.45, size = 2) +
    # geom_point(data = symptom_onset,
    #            inherit.aes = FALSE, aes(x = time, y = ct_value),
    #            colour = custom_pallete[5],
    #            fill = custom_pallete[5],
    #            shape = 2) +
    # annotate(geom = "text", x = 4.5, y = 8, label = 'symptom onset', size = 2) +
    geom_hline(yintercept = 40, linetype = 2, colour = 'gray60') + 
    
    # scale_shape_manual(values = c(16, 17),
    #                    guide = guide_legend(override.aes = list(colour = c(CONFIG$cols[4],
    #                                                                        CONFIG$cols[1]))),
    #                    name = 'Observation',
    #                    labels = c('lab ct', 'symptoms onset')) +
    lims(x = c(-5, 21)) + 
    scale_y_reverse(
      expand = c(0, 0),
      limits = c(45, 15),
      breaks = sort(c(seq(10, 40, 5))),
      labels = c(10,'',20,'', 30 ,'', 40),
      "Ct"
      ) +
    scale_colour_manual(values = custom_pallete[c(3, 4, 7, 8)],
                        labels = c("Covariate 1 Ct result",
                                   "Covariate 2 Ct result"),
                        na.translate = F) +
    labs(x = "Time since first positive test (days)",
         y = "Ct value") +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          legend.position = c(0.9, 0.95),
          legend.background = element_blank(),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key = element_rect(color = NA, fill = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.y = unit(0.01, "cm")) +
    labs(tag = 'A')
  
  return(figure_1_aa)
  
}
