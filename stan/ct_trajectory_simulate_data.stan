#include ct_trajectory_functions.stan

data { 
  int <lower = 0> P;
  real c_0;
  real c_lod;
  real t_e;
  real <lower = 0> lmean;
  real <lower = 0> lsd;
}

parameters {
  // Inferred time of infection
  vector <lower = 0> [P] T_e;
  // real T_e;

  // Timing of p and s
  real <lower = 0> t_p_mean;
  real <lower = 0> t_p_var;
  vector <lower = 0> [P] t_p_raw;

  real <lower = 0> t_s_mean;
  real <lower = 0> t_s_var;
  vector <lower = 0> [P] t_s_raw;

  // Time viral load hits lower limit of detection
  real <lower = 0> t_lod_mean;
  real <lower = 0> t_lod_var;
  vector <lower = 0> [P] t_lod_raw;

  // Ct value of viral load p
  // real <lower = 0, upper = c_lod> c_p_mean;
  real <lower = 0> c_p_var;
  vector <lower = 0> [P] c_p_raw;

  // Ct value at s
  // real <lower = 0, upper = c_lod> c_s_mean;
  real <lower = 0> c_s_var;
  vector <lower = 0> [P] c_s_raw;

  // testing ordered parameters for the two inferred Ct values
  ordered[2] ct_params;

  // Variance parameter for oobservation model
  real <lower = 0> sigma_obs;
}

transformed parameters {
  
  vector [P] t_p;
  vector [P] t_s;
  vector [P] t_lod;
  vector [P] c_s;
  vector [P] c_p;
  vector [P] t_lod_abs;

  // non-centred, hierarchical parameterisation
  t_p = exp(t_p_mean + t_p_var * t_p_raw);
  t_s = exp(t_s_mean + t_s_var * t_s_raw);
  t_lod = exp(t_lod_mean + t_lod_var * t_lod_raw);
  c_s = c_lod * inv_logit(ct_params[1] + c_s_var * c_s_raw);
  c_p = c_s .* inv_logit(ct_params[2] + c_p_var * c_p_raw);
  t_lod_abs = t_p + t_s + t_lod;
}

model {
  // Prior over possible infection times
  T_e ~ cauchy(0, 5);

  // // Viral load peak timing
  t_p_mean ~ cauchy(log(5), 1);
  t_p_var ~ cauchy(0, 1);
  t_p_raw ~ normal(0, 1);

  t_s_mean ~ cauchy(log(10), 2);
  t_s_var ~ cauchy(0, 1);
  t_s_raw ~ normal(0, 1);
  
  // // // Time dropping below limit of detection
  t_lod_mean ~ cauchy(log(15), 2);
  t_lod_var ~ cauchy(0, 1);
  t_lod_raw ~ normal(0, 1);
  
  // // // Ct value at peak
  // c_p_mean ~ normal(0, 1);
  c_p_var ~ cauchy(0, 1);
  c_p_raw ~ normal(0.3, 1);
  
  // // // Ct value at switch to long wane
  // c_s_mean ~ normal(0, 1);
  c_s_var ~ cauchy(0, 1);
  c_s_raw ~ normal(0, 1);

  // testing ordered parameters for the two inferred Ct parameters
  ct_params ~ normal(0, 1);

  // Variation in observation model
  sigma_obs ~ cauchy(0, 5);
}

generated quantities {
  // vector[20] ct;
  // matrix[P, 30] ct;
  // for(i in 1:P) {
  //   for(j in 1:30) {
  //     ct[i, j] = ct_hinge(j, c_0, c_p[i], c_s[i], c_lod, t_e, t_p[i], t_s[i], t_lod_abs[i]);
  //   }
  // }  
  vector[201] ct;
  real k;
  for(j in 1:201) {
    k = (j * 0.1) - 0.1;
    // for(j in 1:20) {
      // ct[j] = ct_hinge(k, c_0, c_p_mean, c_s_mean, c_lod, t_e, t_p_mean, t_s_mean, t_lod_mean);
      ct[j] = ct_hinge_long(k, c_0, ct_params[1], ct_params[2], c_lod, t_e, t_p_mean, t_s_mean, t_lod_mean);
    }
}
