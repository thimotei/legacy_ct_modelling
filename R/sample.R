sample_pop_priors <- function(n_samples,
                              c_lod = 40,
                              switch = FALSE,
                              data_format = "long",
                              scale_type = "natural") {
  
  dt_proc <- data.table(
    c_0 = truncnorm::rtruncnorm(
      n_samples, a = c_lod, mean = c_lod + 10, sd = 2),
    c_p = rnorm(n_samples, 0, 1),
    t_p = rnorm(n_samples, 1.61, 0.5),
    t_lod = rnorm(n_samples, 2.3, 0.5))
  
  if(switch == TRUE) {
    dt_proc[, c_s := rnorm(n_samples, 0, 1)]
    dt_proc[, t_s := rnorm(n_samples, 1.61, 0.5)]
  }
  
  if(scale_type == "natural") {
    dt_proc <- transform_to_natural(dt_proc)
  }
  
  if(data_format == "long") {
    dt_out <- melt(dt_proc,
                   measure.vars = colnames(dt_proc),
                   variable.name = "parameter")
  } else {
    dt_out <- dt_proc
  }
  return(dt_out)
}
