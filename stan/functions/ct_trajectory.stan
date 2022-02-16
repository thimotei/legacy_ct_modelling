functions {

  real ct_hinge_single(real t, real c0, real cp,
                       real clod, real te, real tp,
                       real tlod) {
    real y;
    
    if(t <= te)
    {
      y = c0;
    }
    else if(t > te && t <= te + tp) 
    {
      y = ((t - te)*(cp - c0))/tp  + c0;
    }
    else if(t > te + tp && t <= te + tp + tlod)
    {
      y = ((t - te - tp)*(clod - cp))/tlod  + cp;
    }
    else if(t > tlod)
    {
      y = clod;
    }
    
    return(y);
  }

  real ct_hinge_long(real t, real c0, real cp, real cs, 
                   real clod, real te, real tp, real ts,
                   real tlod)
  {
    real y;
    
    if(t <= te)
    {
      y = c0;
    }
    else if(t > te && t <= te + tp) 
    {
      y = ((t - te)*(cp - c0))/tp  + c0;
    }
    else if(t > te + tp && t <= te + tp + ts)
    {
      y = ((t - te - tp)*(cs - cp))/ts  + cp;
    }
    else if(t > te + tp + ts && t <= te + tp + ts + tlod)
    {
      y = ((t - te - tp - ts)*(clod - cs))/tlod + cs;
    }
    else if(t > tlod) 
    {
      y = clod;  
    }
    
    return(y);
  }
  
  vector ct_hinge_vec_new(vector t, real c0, vector cp, vector cs, real clod, 
                          real te, vector tp, vector ts, vector tlod, int [] patient_ID)
  {
    int N = num_elements(t);
    vector[N] ret;

    for(k in 1:N){
      ret[k] = ct_hinge_long(t[k], c0, 
                             cp[patient_ID[k]], cs[patient_ID[k]],
                             clod, te, tp[patient_ID[k]],
                             ts[patient_ID[k]], tlod[patient_ID[k]]);
    }

    return(ret);
  }
}