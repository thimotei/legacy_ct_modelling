  real piecewise_ct_single(real t, real c0, real cp,
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

  real piecewise_ct(real t, real c0, real cp, real cs, 
                    real clod, real te, real tp, real ts,
                    real tlod) {
    real y;
    real adj_t = t - te;
    if (adj_t <= 0) {
      y = c0;
    }else {
      if (adj_t <= tp) {
        y = (adj_t * (cp - c0)) / tp  + c0;
      }else{
        adj_t = adj_t - tp;
        if (adj_t <= ts) {
          y = (adj_t * (cs - cp)) / ts  + cp;
        }else{
          adj_t = adj_t - ts;
          if (adj_t <= tlod) {
            y = (adj_t * (clod - cs)) / tlod + cs;
          }else{
            y = clod;
          }
        }
      }
    }  
    return(y);
  }
  
  vector piecewise_ct_vec(vector t, real c0, vector cp, vector cs, real clod, 
                         real te, vector tp, vector ts, vector tlod,
                         array[] int patient_ID) {
    int N = num_elements(t);
    vector[N] ret;

    for(k in 1:N){
      ret[k] = piecewise_ct(t[k], c0, 
                            cp[patient_ID[k]], cs[patient_ID[k]],
                            clod, te, tp[patient_ID[k]],
                            ts[patient_ID[k]], tlod[patient_ID[k]]);
    }

    return(ret);
  }
