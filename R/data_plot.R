raw_data_plot <- function(dt){
  
  obs <- data.table::copy(dt)
  
  obs[, combo.id := paste0(id,".",infection_id)]
  obs[, clean.id := .GRP, by = "combo.id"][, clean.id := as.factor(clean.id)]
  obs[, last_swab := max(time_since_first_pos), by = c("clean.id")]
  
  plot_data <- obs[time_since_first_pos >= 0]
  plot_data[, clean.id := forcats::fct_reorder(clean.id, last_swab)]
  
  p <- plot_data %>%
    ggplot(aes(x = time_since_first_pos, 
               y = clean.id,
               col = ct_value)) +
    # geom_line() +
    geom_point(shape= 15, size = 2.5) +
    scale_color_gradient(high = "black", low = "firebrick1", name = "Ct value") +
    scale_x_continuous(breaks = 0:30) +
    cowplot::theme_minimal_grid() +
    # geom_text(aes(y = forcats::fct_reorder(clean.id, last_swab), x = -1 , label = substr(VOC,1,1)), inherit.aes = FALSE, size = 2.5) +
    geom_text(aes(y = forcats::fct_reorder(clean.id, last_swab), x = -0.5, label = paste0("(",no_vaccines,")")), inherit.aes = FALSE, size = 2.5) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
    labs(x = "Days since first positive test") +
    new_scale_color() +
    geom_point(aes(y = clean.id, x = -1, col = VOC), shape = 15)
  
  return(p)
}


min_ct_plot <- function(dt){
  
  obs <- data.table::copy(dt)
  
  obs[, iid := .GRP, c("id", "infection_id")]
  
  pdat <- obs[t >=0 , .(iid, t, age.group, no_vaccines, VOC, symptoms, ct_value, time_since_last_dose, t_onset = t - onset_time)]
  
  pdat <- pdat[, .(no_vaccines = no_vaccines[1], 
                   voc = VOC[1], 
                   symptoms = symptoms[1], 
                   min_ct = min(ct_value), 
                   tsld = time_since_last_dose[1],
                   age.group = age.group[1]), "iid"]
  
  
  p <- ggplot(data = pdat, aes(x = age.group, y = min_ct, fill = symptoms)) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2/3, dotsize = 1.5,
                 position = position_dodge(0.5)) + 
    facet_grid(voc ~ no_vaccines) +
    theme_light() + 
    labs(x = "Age group", y = "Minimum observed Ct value") +
    scale_fill_discrete(name = "")
  
  return(p)
}