plot_obs_ct <- function(ct_dt, ct_traj, pp, traj_alpha = 0.02, onsets = TRUE,
                        clod = 40) {

  if (!missing(pp)) {
    ct_dt <- cbind(
      ct_dt,
      data.table::copy(pp)[, c("t", "id", "pcr_res", "obs") := NULL]
    )
  }

  plot <- ggplot(ct_dt) +
    aes(x = t, y = ct_value, colour = factor(swab_type)) +
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
        aes(ymin = lo90, ymax = hi90, y = NULL), size = 1.1, alpha = 0.2
      ) +
      geom_linerange(
        aes(ymin = lo60, ymax = hi60, y = NULL), size = 1.1, alpha = 0.2
      )
   }

  plot <- plot +
    geom_point(alpha = 0.6)

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
    labs(
      colour = "Swab type", x = "Days since first positive test",
      y = "CT value"
    )
  return(plot)
}

plot_ct_pp <- function(pp, sum_pp, onsets = TRUE, clod = 40, alpha = 0.05,
                       ...) {

  plot <- ggplot(pp) +
    aes(x = t, y = ct_value, group = interaction(.iteration, .chain), ...) +
    scale_colour_brewer(palette = "Dark2")

  if (!is.null(pp$onset_time) & onsets) {
    plot <- plot +
      geom_vline(aes(xintercept = onset_time), linetype = 2, alpha = 0.8)
  }

  if (!is.null(clod)) {
    plot <- plot +
      geom_hline(yintercept = clod, linetype = 3, alpha = 0.8)
  }

   if (!missing(sum_pp)) {
     plot <- plot +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo90, ymax = hi90, y = NULL, group = NULL), alpha = 0.2
      ) +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo60, ymax = hi60, y = NULL, group = NULL), alpha = 0.2
      )
   }

  plot <- plot +
    geom_line(alpha = alpha)

   plot <- plot +
    custom_plot_theme() +
    labs(
      x = "Days since infection",
      y = "CT value"
    )
  return(plot)
}

plot_density <- function(draws, ...) {
  plot <- ggplot(draws) +
    aes(x = value, ...) +
    geom_density(alpha = 0.2) +
    facet_wrap(~variable, nrow = 2, scales = "free") +
    labs(x = "", y = "Probability density")
  return(plot)
}
