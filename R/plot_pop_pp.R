plot_pop_pp <- function(pop_pp_summary_dt) {
  
  pop_pp_out <- pop_pp_summary_dt %>% 
    ggplot() + 
    geom_line(aes(x = time, 
                  y = median, 
                  group = voc,
                  colour = voc),
              linetype = "dashed", alpha = 0.5) +
    geom_line(data = post_pred_samples[iteration %in% 1:100],
              aes(x = time, 
                  y = value, 
                  group = interaction(voc, iteration), 
                  colour = voc),
              alpha = 0.05) +
    geom_ribbon(aes(x = time,
                    ymin = lo90,
                    ymax = hi90, 
                    fill = voc),
                alpha = 0.1) +
    theme(legend.position = "none") + 
    labs(x = "Time (days)", y = "Ct value")
  
  return(pop_pp_out)
  
}
