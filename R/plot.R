plot_obs_ct <- function(ct_dt, ct_traj, pp, traj_alpha = 0.02, onsets = TRUE,
                        clod = 40) {
  plot <- ggplot(ct_dt) +
    aes(x = t, y = ct_value, colour = factor(pcr_res)) +
    scale_colour_brewer(palette = "Dark2")

  if (!is.null(ct_dt$onset_time) & onsets) {
    plot <- plot +
      geom_vline(aes(xintercept = onset_time), linetype = 2, alpha = 0.8)
  }

  if (!is.null(clod)) {
    plot <- plot +
      geom_hline(yintercept = clod, linetype = 3, alpha = 0.8)
  }

   if (!missing(pp)) {
     plot <- plot +
      geom_linerange(
        data = pp, aes(ymin = lo, ymax = hi, y = NULL), size = 1.1, alpha = 0.2
      )
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
