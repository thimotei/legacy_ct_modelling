ct_hinge_dt <- function(c_peak, t_peak, c_switch, t_switch, t_LOD){
  
  out <- vapply(X = seq(0, 40, 0.1), FUN = ct_hinge, FUN.VALUE = 1, c_zero = c_zero,
                c_switch = c_switch, t_switch = t_switch, t_eclipse = 0,
                c_LOD = c_lod, c_peak = c_peak, t_peak = t_peak, t_LOD = t_lod)
  
  # dtout <- data.table(value = out)
  # dtout$vocl <- rep(voc, nrow(dtout))
  # dtout$tsi <- seq(0, 40, 0.1)
  
  return(out)
}

ct_hinge <- function(a, c_zero, 
                     c_peak, c_switch, 
                     c_LOD, t_eclipse, 
                     t_peak, t_switch, t_LOD) {
  
  if(a <= t_eclipse){
    out = c_zero
  }else if(a > t_eclipse){
    if(a <= (t_eclipse + t_peak)){
      out = c_zero + ((c_peak - c_zero) / t_peak) * (a - t_eclipse)
    }else if(a > (t_eclipse + t_peak)){
      if(a <= (t_eclipse + t_peak + t_switch)){
        out = c_peak + ((c_switch - c_peak) / t_switch) * (a - t_eclipse - t_peak)
      }else if(a > (t_eclipse + t_peak + t_switch)){
        out = c_switch + ((c_LOD - c_switch) / (t_LOD - t_switch - t_peak - t_eclipse)) * (a - t_eclipse - t_peak - t_switch)
      }
    }
  }
  
  return(out)
}