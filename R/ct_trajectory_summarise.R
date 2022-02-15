ct_trajectory_summarise <- function(draws, by = c("id", "time")) {

out <- draws[,
    .(me = quantile(value, c(0.5)), lo = quantile(value, c(0.025)), 
      hi = quantile(value, c(0.975))
    ),
    by = by
  ]
  return(out)
}
