#include ct_trajectory_functions.stan

data {
  int P; // number of patients
  int N; // number of tests
  real c_0; // Ct value before detection
  real c_lod; // Ct value at limit of detection 
  real t_e; 
  real lmean; // mean of incubation period used
  real lsd; // standard deviation of incubation period used
  int id[N]; // id of person
  int pcr_res[N]; // boolean test result
  vector[N] day_rel; // day of test (integer)
  vector[N] ct_value; // Ct value of test
  // vector [P] symp_rel;
  // vector [P] te_upper_bound; // upper bound on infection time
}

parameters {
  // Inferred time of infection
  vector[P] T_e;
  
  // Hyperparameters
  // Ct value of viral load p
  real c_p_mean;
  real <lower = 0> c_p_var;
  vector [P] c_p_raw;

  // Ct value at s
  real c_s_mean;
  real <lower = 0> c_s_var;
  vector [P] c_s_raw;
  
  // Timing of peak
  real t_p_mean;
  real<lower=0> t_p_var;
  vector[P] t_p_raw;

  // Timing of switch
  real t_s_mean;
  real <lower = 0> t_s_var;
  vector [P] t_s_raw;

  // Time viral load hits lower limit of detection
  real t_lod_mean;
  real <lower = 0> t_lod_var;
  vector [P] t_lod_raw;

  // Variance parameter for oobservation model
  real<lower=0> sigma_obs;
}

transformed parameters {

  vector[P] t_p;
  vector[P] t_s;
  vector[P] t_lod;
  vector[P] c_p;
  vector[P] c_s;
  vector[P] t_lod_abs;

  // individual-level parameters
  // non-centred, hierarchical parameterisation
  t_p = exp(t_p_mean + t_p_var * t_p_raw);
  t_s = exp(t_s_mean + t_s_var * t_s_raw);
  t_lod = exp(t_lod_mean + t_lod_var * t_lod_raw);
  // Parameterise c_switch as proportion of c_LOD
  c_s = c_lod * inv_logit(c_s_mean + c_s_var * c_s_raw);
  // Parameterise c_peak as proportion of c_switch
  c_p = c_s .* inv_logit(c_p_mean + c_p_var * c_p_raw);
  t_lod_abs = t_p + t_s + t_lod;
  
}

model {
  vector[N] diff = day_rel + T_e[id];
  vector[N] exp_ct;

  // Prior over possible infection times
  T_e[i] ~ normal(3, 2);
  
  // Ct value at peak
  c_p_mean ~ normal(0, 1); #mean at 50% of switch value
  c_p_var ~ normal(0, 0.1) T[0,];
  c_p_raw ~ std_normal();

  // Ct value at switch to long wane
  c_s_mean ~ normal(0, 1); # mean at 50% of maximum ct
  c_s_var ~ normal(0, 0.1) T[0,];
  c_s_raw ~ std_normal();

  // Viral load peak timing
  t_p_mean ~ normal(1.61, 0.5); #mean at log(5)
  t_p_var ~ normal(0, 0.1) T[0,];
  t_p_raw ~ normal(0, 1);

  t_s_mean ~ normal(1.61, 0.5); #mean at log(5) + peak timing
  t_s_var ~ normal(0, 0.1) T[0,];
  t_s_raw ~ std_normal();

  // Time dropping below limit of detection
  t_lod_mean ~ normal(2.3, 0.5); # mean at log(10) + peak + scale timing 
  t_lod_var ~ normal(0, 0.1) T[0,];
  t_lod_raw ~ std_normal();

  // // Variation in observation model (% scale of C_lod)
  sigma_obs ~ normal(0, 0.025) T[0,];

  // // component of likelihood for time of exposure
  // for(j in 1:P) {
  //   // likelihood for time of exposure using the CDF of the incubation period
  //   // and known symptom onset times. Need to form a tiny artificial
  //   // time-window to centre the estimates properly (using 0.05 of a day).
  //   // Is there a better/less hacky way than this?
  //   target += log(
  //     lognormal_cdf(symp_rel[j] - T_e[j], lmean, lsd) -
  //     lognormal_cdf(symp_rel[j] - 1 - T_e[j], lmean, lsd));
  // }

  // Expected ct value given viral load parameters
  exp_ct = ct_hinge_vec_new(diff, c_0, c_p, c_s, c_lod, t_e, t_p, t_s, 
                            t_lod_abs, id);

  // component of likelihood for expected ct values
  for(j in 1:N) {
    // If positive result: P(observed ct | expected ct)*P(Ct detected | expected ct)
    if(pcr_res[j] == 1) {
      ct_value[j] ~ normal(exp_ct[j], sigma_obs) T[0, c_lod];
      target += normal_lcdf(c_lod | exp_ct[j], sigma_obs);
      }
    // if negative result: P(Ct not detected | expected ct)
    else if(pcr_res[j] == 0) {
      target += normal_lccdf(c_lod | exp_ct[j], sigma_obs);
    }
  }

}

generated quantities {
  matrix[P, 61] ct;
  for(i in 1:P) {
    for(j in 1:61) {
      ct[i, j] = ct_hinge_long(j - 1, c_0, c_p[i], c_s[i], c_lod, t_e, t_p[i], t_s[i], t_lod_abs[i]);
    }
  }
}
