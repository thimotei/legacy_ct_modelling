#include ct_trajectory_functions_old.stan

data {
  int <lower = 0> N; // number of tests
  int <lower = 0> P; // number of patients
  int id[N]; // id of person
  vector [N] day_rel; // day of test (integer)
  real ct_value[N]; // Ct value of test
  int pcr_res[N]; // boolean test result
  // vector [P] symp_rel;
  // vector [P] te_upper_bound; // upper bound on infection time
  real c_0;
  real c_lod;
  real t_e;
  real <lower = 0> lmean;
  real <lower = 0> lsd;
}

parameters {
  // Inferred time of infection
  vector <lower = 0> [P] T_e;

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
  real <lower = 0> c_p_mean;
  real <lower = 0> c_p_var;
  vector <lower = 0> [P] c_p_raw;

  // Ct value at s
  real <lower = 0> c_s_mean;
  real <lower = 0> c_s_var;
  vector <lower = 0> [P] c_s_raw;
  
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
  c_s = c_lod * inv_logit(c_s_mean + c_s_var * c_s_raw);
  c_p = c_s .* inv_logit(c_p_mean + c_p_var * c_p_raw);
  t_lod_abs  = t_p + t_s + t_lod;

}

model {
  vector [N] diff = day_rel - T_e[id];
  // vector [N] diff = day_rel - T_e;
  vector [N] exp_ct;

  // component of likelihood for time of exposure
  // for(j in 1:P) {
  //   // likelihood for time of exposure using the CDF of the incubation period
  //   // and known symptom onset times. Need to form a tiny artificial
  //   // time-window to centre the estimates properly (using 0.05 of a day).
  //   // Is there a better/less hacky way than this?
  //   target += log(
  //     lognormal_cdf(symp_rel[j] - T_e[j], lmean, lsd) - 
  //     lognormal_cdf(symp_rel[j] - 0.05 - T_e[j], lmean, lsd));
  // }
  
  // Expected ct value given viral load parameters
  exp_ct = ct_hinge_vec(diff, c_0, c_p, c_s, c_lod, t_e, t_p, t_s, t_lod_abs, id);

  // component of likelihood for expected ct values
  for(j in 1:N) {
    // If positive result: P(observed ct | expected ct)*P(Ct detected | expected ct)
    if(pcr_res[j] == 1) {
      ct_value[j] ~ normal(exp_ct[j], sigma_obs) T[, c_lod];
      target +=  normal_lcdf(c_lod | exp_ct[j], sigma_obs);
      }
    // if negative result: P(Ct not detected | expected ct)
    else {
      target += normal_lccdf(c_lod | exp_ct[j], sigma_obs);
    }
  }

  // Prior over possible infection times
  T_e ~ cauchy(0, 1);

  // Viral load peak timing
  t_p_mean ~ cauchy(log(5), 1);
  t_p_var ~ cauchy(0, 1);
  t_p_raw ~ normal(0, 1);

  t_s_mean ~ cauchy(log(8), 1);
  t_s_var ~ cauchy(0, 1);
  t_s_raw ~ normal(0, 1);

  // // Time dropping below limit of detection
  t_lod_mean ~ cauchy(log(12), 1);
  t_lod_var ~ cauchy(0, 1);
  t_lod_raw ~ normal(0, 1);

  // // Ct value at peak
  c_p_mean ~ cauchy(log(0.1), 1);
  c_p_var ~ cauchy(0, 1);
  c_p_raw ~ normal(0.3, 1);

  // // Ct value at switch to long wane
  c_s_mean ~ cauchy(log(0.3), 1);
  c_s_var ~ cauchy(0, 1);
  c_s_raw ~ normal(0, 1);

  // // Variation in observation model
  sigma_obs ~ cauchy(0, 2);
}

generated quantities {
  matrix[P, 30] ct;
  for(i in 1:P) {
    for(j in 1:30) {
      ct[i, j] = ct_hinge(j, c_0, c_p[i], c_s[i], c_lod, t_e, t_p[i], t_s[i], t_lod_abs[i]);
    }
  }
}