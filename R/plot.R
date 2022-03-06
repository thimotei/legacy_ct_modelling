custom_plot_theme <- function(flip = FALSE, legend_arg = FALSE) {

  custom_plot_theme <- list(
    theme_bw(11),
    theme(strip.placement = "outside"),
    geom_vline(aes(xintercept = -Inf)),
    guides(color = guide_legend(override.aes = list(fill = NA)))
  )

  if (!flip) {
    append(custom_plot_theme, geom_hline(aes(yintercept = -Inf)))
  }

  if (!legend_arg) {
    append(custom_plot_theme,
           theme(legend.title = element_blank()))
  }

  return(custom_plot_theme)
}


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
    aes(x = t, y = ct_value, group = interaction(.iteration, .chain), ...)

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
        aes(ymin = lo90, ymax = hi90, y = NULL, group = NULL, col = NULL), alpha = 0.15
      ) +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo60, ymax = hi60, y = NULL, group = NULL, col = NULL), alpha = 0.15
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

plot_ip_pp <- function(pp, sum_pp, onsets = TRUE, alpha = 0.05, ...) {
  pp <- data.table::copy(pp)[, ct_value := value]
  plot <- plot_ct_pp(
    pp, sum_pp, onsets = TRUE, alpha = alpha, clod = NULL, ...
  ) +
    labs(y = "Probability of symptom onset")
    return(plot)
}

plot_pp_from_fit <- function(fit, obs, samples = 10, alpha = 0.05) {
  ct_draws <- extract_ct_trajectories(fit)

  ct_summary <- summarise_draws(
    data.table::copy(ct_draws)[,
      time_since_first_pos := as.integer(time_since_first_pos)
      ],
    by = c("id", "time_since_first_pos")
  )

  ct_pp <- extract_posterior_predictions(fit, obs)
  ct_pp <- summarise_draws(
    ct_pp[, value := sim_ct], by = c("id", "t", "pcr_res", "obs")
  )

  # Plotting summaries of fitted trajectories against simulated data
  plot <- plot_obs_ct(
    obs, ct_draws[iteration <= ceiling(samples / max(chain))],
    ct_pp, traj_alpha = alpha
  )
  return(plot)
}

plot_density <- function(draws, ...) {
  plot <- ggplot(draws) +
    geom_density(aes(x = value, y = ..scaled.., ...), alpha = 0.2) +
    custom_plot_theme() +
    labs(x = "", y = "Density")
  return(plot)
}

plot_ct_summary <- function(draws, time_range = seq(0, 60, by = 0.25),
                            samples = 100, by = c(), traj_alpha = 0.05,
                            simulated_samples = 1000, ...) {
  pop_draws <- extract_pop_params(draws, by = by)

  pop_ct_draws <- pop_draws[.draw <= simulated_samples] %>%
    transform_to_model() %>%
    simulate_cts(time_range = time_range, obs_noise = FALSE)

  sum_cols <- c("value", "t", by)
  pop_ct_sum <- summarise_draws(
    pop_ct_draws[, value := ct_value][, ..sum_cols],
    by = setdiff(sum_cols, "value")
  )

  ct_pp_plot <- plot_ct_pp(
    pop_ct_draws[.draw <= samples], pop_ct_sum, alpha = traj_alpha, ...
  )

  param_pp_plot <- pop_draws %>%
    transform_to_natural() %>%
    melt_draws(ids = c(".chain", ".iteration", ".draw", by)) %>%
    update_variable_labels() %>%
    plot_density(...) +
    ggplot2::facet_wrap(~variable, nrow = 2, scales = "free_x")

  plot <- param_pp_plot / ct_pp_plot +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  return(plot)
}

plot_ip_summary <- function(draws, time_range = seq(0, 20, by = 0.25),
                            samples = 100, by = c(), traj_alpha = 0.05, simulated_samples = 1000, ...) {
  ip_draws <- extract_ip_params(draws, by = by)

  pop_ip_draws <- ip_draws[.draw <= simulated_samples] %>%
    simulate_ips(time_range = time_range)

  sum_cols <- c("value", "t", by)
  pop_ip_sum <- summarise_draws(
    pop_ip_draws[, ..sum_cols],
    by = setdiff(sum_cols, "value")
  )

  ip_pp_plot <- plot_ip_pp(
    pop_ip_draws[.draw <= samples], pop_ip_sum, alpha = traj_alpha, ...
  )

  param_pp_plot <- ip_draws %>%
    transform_ip_to_natural() %>%
    melt_draws(ids = c(".chain", ".iteration", ".draw", by)) %>%
    update_variable_labels() %>%
    plot_density(...) +
    ggplot2::facet_wrap(~variable, nrow = 2, scales = "free_x")

  plot <- param_pp_plot / ip_pp_plot +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  return(plot)
}

plot_summary <- function(draws, ct_time_range = seq(0, 60, by = 0.25),     
                         ip_time_range = seq(0, 20, by = 0.25), samples = 100,
                         by = c(), traj_alpha = 0.05, simulated_samples = 1000, ...) {
  ct_pp <- plot_ct_summary(
    draws, time_range = ct_time_range, samples = samples, by = by, simulated_samples = simulated_samples, traj_alpha = traj_alpha,
    ...
  )

  ip_pp <- plot_ip_summary(
    draws, time_range = ip_time_range, samples = samples, by = by, simulated_samples = simulated_samples, traj_alpha = traj_alpha,
    ...
  )

  parameter_pp <- ((ct_pp) | (ip_pp)) +
    patchwork::plot_layout(
      guides = "collect", widths = c(3, 2),
    ) &
    theme(legend.position = "bottom")
  return(parameter_pp)
}

plot_effects <- function(effects,  position = "identity", ...) {
 
  eff_plot <- ggplot(effects) +
    aes(y = variable, ...) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_linerange(
      aes(xmin = lo90, xmax = hi90), 
      size = 3, alpha = 0.3,
      position = position
    ) +
    geom_linerange(
      aes(xmin = lo60, xmax = hi60), 
      size = 3, alpha = 0.3,
      position = position
    ) +
    custom_plot_theme() +
    theme(legend.position = "bottom") +
    scale_x_continuous(trans = "log") +
    labs(x = "Effect size", y = "Variable modified")
  return(eff_plot)
}