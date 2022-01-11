functions {
  real ct_hinge(real a, real c_zero, 
  real c_peak, real c_switch, 
  real c_LOD, real t_eclipse, 
  real t_peak, real t_switch, real t_LOD) {
    real out;
    
    if(a <= t_eclipse){
      out = c_zero;
    }else if(a > t_eclipse){
      if(a <= (t_eclipse + t_peak)){
        out = c_zero + ((c_peak - c_zero) / t_peak) * (a - t_eclipse);
      }else if(a > (t_eclipse + t_peak)){
        if(a <= (t_eclipse + t_peak + t_switch)){
          out = c_peak + ((c_switch - c_peak) / t_switch) * (a - t_eclipse - t_peak);
        }else if(a > (t_eclipse + t_peak + t_switch)){
            out = c_switch + ((c_LOD - c_switch) / (t_LOD - t_switch - t_peak - t_eclipse)) * (a - t_eclipse - t_peak - t_switch);
        }
      }
    }
    
    return(out);
  }
  
  vector ct_hinge_vec(vector a, real c_zero, vector c_peak, vector c_switch, real c_LOD, real t_eclipse, vector t_peak, vector t_switch, vector t_LOD, int [] patient_ID) {
    
    int N = num_elements(a);
    vector[N] ret;
    
    for(k in 1:N){
      ret[k] = ct_hinge(a[k], c_zero, c_peak[patient_ID[k]], c_switch[patient_ID[k]], c_LOD, t_eclipse, t_peak[patient_ID[k]], t_switch[patient_ID[k]], t_LOD[patient_ID[k]]);
    }
    
    return(ret);
  }
}