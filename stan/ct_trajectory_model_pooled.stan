#include ct_trajectory_functions.stan

data {
  int <lower = 0> N; // number of tests
  int <lower = 0> P; // number of patients
  int id[N]; // id of person
  vector [N] day_rel; // day of test (integer)
  real ct_value[N]; // Ct value of test
  int pcr_res[N]; // boolean test result
  real c_0;
  real c_lod;
  real t_e;
  real <lower = 0> lmean;
  real <lower = 0> lsd;
}

parameters {
  
  // Inferred time of infection
  real T_e;

  // Timing of peak
  real <lower = 0> t_p_mean;

  // Timing of switch
  real <lower = 0> t_s_mean;

  // Time viral load hits lower limit of detection
  real <lower = 0> t_lod_mean;

  // Ct value of viral load p
  real <lower = 0, upper = c_lod> c_p_mean;

  // Ct value at s
  real <lower = 0, upper = c_lod> c_s_mean;
  
  // Variance parameter for oobservation model
  real <lower = 0> sigma_obs;
}

transformed parameters {
  real  t_lod_abs;
  t_lod_abs = t_p_mean + t_s_mean + t_lod_mean;
}

model {
  
  vector [N] diff = day_rel - T_e;
  vector [N] exp_ct;
  
  // likelihood for time of exposure using the CDF of the incubation period
  // and known symptom onset times
  // for(j in 1:P) {
  //   target += log(
  //     lognormal_cdf(symp_rel[j] - T_e[j], lmean, lsd) -
  //     lognormal_cdf(symp_rel[j] - 1.0 - T_e[j], lmean, lsd));
  // }
  // 
  // Expected ct value given viral load parameters
  for(i in 1:N) {
    exp_ct[i] = ct_hinge_long(diff[i], c_0, c_p_mean, c_s_mean, c_lod, t_e, t_p_mean, t_s_mean, t_lod_abs);
  }

  // component of likelihood for expected ct values
  for(j in 1:N) {
    // If positive result: P(observed ct | expected ct)*P(Ct detected | expected ct)
    if(pcr_res[j] == 1) {
      ct_value[j] ~ normal(exp_ct[j], sigma_obs) T[, c_lod];
      target +=  normal_lcdf(c_lod | exp_ct[j], sigma_obs);
      }
    // if negative result: P(Ct not detected | expected ct)
    else if(pcr_res[j] == 0) {
      target += normal_lccdf(c_lod | exp_ct[j], sigma_obs);
    }
  }

  // latent infection time prior
  // T_e ~ cauchy(log(5), 5);
  // 
  // // semi-mechanistic Ct trajectory model priors
  // t_p_mean ~ cauchy(log(10), 5);
  // t_s_mean ~ cauchy(log(10), 5);
  // t_lod_mean ~ cauchy(log(20), 5);
  // c_p_mean ~ cauchy(0.4, 5);
  // c_s_mean ~ cauchy(0.6, 5);
  
  // latent infection time prior
  T_e ~ normal(0, 5);

  // semi-mechanistic Ct trajectory model priors
  t_p_mean ~ normal(5, 2);
  t_s_mean ~ normal(15, 2);
  t_lod_mean ~ normal(20, 2);
  c_p_mean ~ normal(0.1, 2);
  c_s_mean ~ normal(0.4, 2);
  
  // observation error
  sigma_obs ~ normal(0, 2);
}

generated quantities {

  vector[601] ct;
  real k;
  
  for(j in 1:601) {
    k = (j * 0.1) - 0.1;
      ct[j] = ct_hinge_long(k, c_0, c_p_mean, c_s_mean, c_lod, t_e, t_p_mean, t_s_mean, t_lod_abs);
    }
}
