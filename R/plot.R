plot_obs_ct <- function(ct_dt, ct_traj, traj_alpha = 0.01, onsets = TRUE) {
  plot <- ggplot(ct_dt) +
    aes(x = t, y = ct_value, colour = factor(pcr_res))

  if (!is.null(ct_dt$onset_time) & onsets) {
    plot <- plot +
      geom_vline(aes(xintercept = onset_time), linetype = 2, alpha = 0.8)
  }

  plot <- plot +
    geom_point()

   if (!missing(ct_traj)) {
     plot <- plot +
      geom_line(
        data = ct_traj,
        aes(y = value, x = time_since_first_pos,
            group = interaction(iteration, chain)
        ),
        colour = "black", alpha = traj_alpha
      )
   }

   plot <- plot +
      facet_wrap(vars(factor(id))) +
      custom_plot_theme() +
      theme(legend.position = "bottom") +
      labs(colour = "PCR result", x = "Days since first positive test",
           y = "CT value")
  return(plot)
}
