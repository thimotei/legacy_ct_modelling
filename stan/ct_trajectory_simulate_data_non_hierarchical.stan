#include ct_trajectory_functions.stan

data { 
  int <lower = 0> P; // number of patients
  int t_max;
  real c_0;
  real c_lod;
  // real t_e;
  real <lower = 0> lmean;
  real <lower = 0> lsd;
}

parameters {
  // Inferred time of infection
  // real <lower = 0> T_e;
  real <lower = 0> t_e;
  real <lower = 0> t_p;
  real <lower = 0> t_s;
  real <lower = 0> t_lod;
  real <lower = 0, upper = c_lod> c_p;
  real <lower = 0, upper = c_lod> c_s;
  real <lower = 0> sigma_obs;
}

model {
  // Prior over possible infection times
  // only running chains for one iteration,
  // so prior values not being used. Stan uses a
  // Uniform (-2, 2) distribution for the first iteration
  // in the absence of user-defined initial conditions
}

generated quantities {
  
  vector[t_max] ct;
  
  for(t in 1:t_max) {
      ct[t] = ct_hinge_long(t - 1, c_0, c_p, c_s, c_lod, t_e, t_p, t_s, t_lod);
    }
    
}
