ct_hinge_long <- function(t, c0, cp, cs, clod, te, tp, ts, tlod) {
  
  if(t <= te)
  {
    y = c0
  }
  else if(t > te && t <= te + tp) 
  {
    y = ((t - te)*(cp - c0))/tp  + c0
  }
  else if(t > te + tp && t <= te + tp + ts)
  {
    y = ((t - te - tp)*(cs - cp))/ts  + cp
  }
  else if(t > te + tp + ts && t <= te + tp + ts + tlod)
  {
    y = ((t - te - tp - ts)*(clod - cs))/tlod + cs
  }
  else if(t > tlod) 
  {
    y = clod
  }
  
  return(y)
}