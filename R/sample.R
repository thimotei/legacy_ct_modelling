sample_from_pop_priors <- function(samples, 
                                   c_lod = 40,
                                   predictor_arg) {
  
  dt_wide <- data.table(
    c_0 = rtruncnorm(samples, a = c_lod, c_lod + 10, 5),
    c_p = rnorm(samples, 0, 1),
    t_p = rnorm(samples, 1.61, 0.5),
    t_lod = rnorm(samples, 2.3, 0.5)
  )
  
  dt_natural_units <- transform_to_model(dt_wide)
  
  dt_long <- melt(
    dt_natural_units, 
    measure.vars = c("c_0", "c_p", "c_s",
                     "t_p",  "t_s", "t_lod"),
    variable.name = "parameter")[,
                                 type := "Prior"]
  
  dt_out <- dt_long[, predictor := predictor_arg]
  
  return(dt_long)
  
}
