plot_obs_ct <- function(ct_dt, ct_traj) {
  plot <- ggplot(ct_dt) +
    aes(x = t, y = ct_value, colour = factor(pcr_res)) +

   if (!missing(ct_traj)) {
     plot <- plot +
      geom_line(data = ct_traj, aes(y = me, x = time_since_first_pos)) +
      geom_ribbon(
        data = ct_traj,
        aes(ymin = lo, ymax = hi, x = time_since_first_pos, colour = NA),
        alpha = 0.2
      )
   }

   plot <- plot +
      facet_wrap(vars(factor(id))) +
      custom_plot_theme()
  return(plot)
}
