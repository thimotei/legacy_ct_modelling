summarise_draws <- function(draws, by = c("id", "time")) {

out <- draws[,
    .(median = quantile(value, c(0.5)),
      lo90 = quantile(value, c(0.05)),
      lo60 = quantile(value, c(0.20)),
      hi60 = quantile(value, c(0.80)),
      hi90 = quantile(value, c(0.95))
    ),
    by = by
  ]
  return(out)
}
