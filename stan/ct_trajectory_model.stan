functions{ 
#include functions/ct_trajectory.stan
#include functions/combine_effects.stan
#include functions/truncated_normal_rng.stan
#include functions/censor.stan
}

data {
  int P; // number of patients
  int N; // number of tests
  real c_lod; // Ct value at limit of detection 
  real t_e; 
  real lmean[2]; // mean of incubation period used (+ sd)
  real lsd[2]; // standard deviation of incubation period used (+ sd)
  int id[N]; // id of person
  int pcr_res[N]; // boolean test result
  vector[N] day_rel; // day of test (integer)
  vector[N] ct_value; // Ct value of test
  int swab_types; // Number of swab types used
  int swab_type[N]; // Swab type per sample
  int any_onsets;
  vector[P] onset_avail;
  vector[P] onset_time;
  int likelihood;
  int preds; // Number of predictors
  real preds_sd; // Standard deviation of predictor coeffs
  matrix[P, preds + 1] design; //Design matrix
  int adj_t_p; // Should time at peak be adjusted
  int adj_t_s; // Should time at switch be adjusted
  int adj_t_lod; // Should time at LOD be adjusted
  int adj_c_p; // Should CT at peak be adjusted
  int adj_c_s; // Should CT at switch be adjusted
  int adj_inc_mean; // Should incubation period mean be adjusted
}

transformed data {
  vector[P] T_e_bound;
  for (i in 1:P) {
    T_e_bound[i] = max({-onset_time[i], 0});
  }
}

parameters {
  // Inferred time of infection
  vector<lower = T_e_bound>[P] T_e;
  
  //Incubation period
  real inc_mean[any_onsets];
  real<lower = 0> inc_sd[any_onsets];
  
  // Ct value before detection
  real<lower = c_lod> c_0; 

  // Hyperparameters
  // Ct value of viral load p
  real c_p_mean;
  real<lower = 0>c_p_var;
  vector[P] c_p_raw;

  // Ct value at s
  real c_s_mean;
  real<lower = 0> c_s_var;
  vector[P] c_s_raw;
  
  // Timing of peak
  real t_p_mean;
  real<lower = 0> t_p_var;
  vector[P] t_p_raw;

  // Timing of switch
  real t_s_mean;
  real<lower = 0> t_s_var;
  vector[P] t_s_raw;

  // Time viral load hits lower limit of detection
  real t_lod_mean;
  real<lower = 0> t_lod_var;
  vector[P] t_lod_raw;
  
  // Swab type intercept and gradient
  vector[swab_types] swab_type_int;
  vector[swab_types] swab_type_grad;

  // Variance parameter for oobservation model
  real<lower = 0> sigma;

  // Coefficients
  vector[preds && adj_t_p ? preds : 0] beta_t_p;
  vector[preds && adj_t_s ? preds : 0] beta_t_s;
  vector[preds && adj_t_lod ? preds : 0] beta_t_lod;
  vector[preds && adj_c_s ? preds : 0] beta_c_s;
  vector[preds && adj_c_p ? preds : 0] beta_c_p;
  vector[preds && adj_inc_mean ? preds : 0] beta_inc_mean;
}

transformed parameters {
  vector[P] t_p;
  vector[P] t_s;
  vector[P] t_lod;
  vector[P] c_p;
  vector[P] c_s;
  vector[P] t_lod_abs;
  vector[N] diff;
  vector[N] exp_ct;
  vector[swab_types + 1] st_int;
  vector[swab_types + 1] st_grad;
  vector[N] adj_exp_ct;
  // individual-level time since infection parameters
  t_p = combine_effects(t_p_mean, beta_t_p, design);
  t_p = exp(t_p + t_p_var * t_p_raw);

  t_s = combine_effects(t_s_mean, beta_t_s, design);
  t_s = exp(t_s + t_s_var * t_s_raw);

  t_lod = combine_effects(t_lod_mean, beta_t_lod, design);
  t_lod = exp(t_lod + t_lod_var * t_lod_raw);

  // Parameterise c_switch as proportion of c_0
  c_s = combine_effects(c_s_mean, beta_c_s, design);
  c_s = c_0 * inv_logit(c_s + c_s_var * c_s_raw);

  // Parameterise c_peak as proportion of c_switch
  c_p = combine_effects(c_p_mean, beta_c_p, design);
  c_p = c_s .* inv_logit(c_p + c_p_var * c_p_raw);

  // Make times absolute
  t_lod_abs = t_p + t_s + t_lod;
  // Adjust observed times since first test to be time since infection
  diff = day_rel + T_e[id];

  // Expected ct value given viral load parameters
  exp_ct = ct_hinge_vec_new(diff, c_0, c_p, c_s, c_0, t_e, t_p, t_s, 
                            t_lod_abs, id);

  // Adjust Swab types
  st_int[1] = 0;
  st_grad[1] = 1;
  if (swab_types) {
    st_int[2:(swab_types + 1)] = swab_type_int;
    st_grad[2:(swab_types + 1)] = swab_type_grad;
  }
  adj_exp_ct = st_int[swab_type] + st_grad[swab_type] .* exp_ct;
}

model {
  // Prior over possible infection times relative to first
  // positive test or symtom onset.
  // Assumes that the first positive test is not a false positive.
  for (i in 1:P) {
    T_e[i] ~ normal(T_e_bound[i] + 5, 5) T[T_e_bound[i],];
  }
  // CT value prior/post detection
  c_0 ~ normal(c_lod + 10, 5) T[0, ];
  
  // Ct value at peak
  c_p_mean ~ normal(0, 1); //mean at 50% of switch value
  c_p_var ~ normal(0, 0.25) T[0,];
  c_p_raw ~ std_normal();

  // Ct value at switch to long wane
  c_s_mean ~ normal(0, 1); //mean at 50% of maximum ct
  c_s_var ~ normal(0, 0.25) T[0,];
  c_s_raw ~ std_normal();

  // Viral load peak timing
  t_p_mean ~ normal(1.61, 0.5); //mean at log(5)
  t_p_var ~ normal(0, 0.25) T[0,];
  t_p_raw ~ std_normal();

  t_s_mean ~ normal(1.61, 0.5); //mean at log(5) + peak timing
  t_s_var ~ normal(0, 0.25) T[0,];
  t_s_raw ~ std_normal();

  // Time dropping below limit of detection
  t_lod_mean ~ normal(2.3, 0.5); //mean at log(10) + peak + scale timing 
  t_lod_var ~ normal(0, 0.25) T[0,];
  t_lod_raw ~ std_normal();

  // If multiple swab types make linear adjustments
  if (swab_types) {
    swab_type_int ~ std_normal();
    swab_type_grad ~ normal(1, 1);
  }

  // Variation in observation model
  sigma ~ normal(0, 2) T[0,];

  // Coefficients priors for predictors
  if (preds) {
    if (adj_t_p) {
      beta_t_p ~ normal(0, preds_sd);
    }
    if (adj_t_s) {
      beta_t_s ~ normal(0, preds_sd);
    }
    if (adj_t_lod) {
      beta_t_p ~ normal(0, preds_sd);
    }
    if (adj_c_p) {
      beta_c_p ~ normal(0, preds_sd);
    }
    if (adj_c_s) {
      beta_c_s ~ normal(0, preds_sd);
    }
    if (adj_inc_mean) {
      beta_inc_mean ~ normal(0, preds_sd);
    } 
  }

  if (any_onsets && likelihood) {
    vector[P] inc_mean_p;

    // Priors on the incubation period
    inc_mean[1] ~ normal(lmean[1], lmean[2]);
    inc_mean_p = combine_effects(inc_mean[1], beta_inc_mean, design);
    inc_sd[1] ~ normal(lsd[1], lsd[2]) T[0, ];
   // component of likelihood for time of exposure
   for(j in 1:P) {
   // likelihood for time of exposure using the CDF of the incubation period
   // and known symptom onset day
   // What is the probability onsets on observed day
    if (onset_avail[j]) {
      real onset_from_inf = onset_time[j] + T_e[j];
      real onset_window = max({0, onset_from_inf - 1});
        target += log_diff_exp(
          lognormal_lcdf(onset_from_inf | inc_mean_p[j], inc_sd[1]),
          lognormal_lcdf(onset_window | inc_mean_p[j], inc_sd[1])
        );
     }
    }
  }

  if (likelihood) {
    // component of likelihood for expected ct values
    for(j in 1:N) {
      // If non-censored result: P(observed ct | expected ct)
      // Truncated above 0 and below latent limit of detection
      if(pcr_res[j]) {
        ct_value[j] ~ normal(adj_exp_ct[j], sigma) T[0, c_0];
      } else{
      // if censored result: P(Ct not detected | expected ct)
        target += normal_lccdf(c_lod | adj_exp_ct[j], sigma);
      }
    }
  }
}

generated quantities {
  matrix[P, 61] ct;
  vector[N] sim_ct;
  for (i in 1:N) {
    sim_ct[i] = truncated_normal_rng(adj_exp_ct[i], sigma, 0, c_0);
    sim_ct[i] = censor(sim_ct[i], c_lod);
  }
  for(i in 1:P) {
    for(j in 1:61) {
      ct[i, j] = ct_hinge_long(j - 1, c_0, c_p[i], c_s[i], c_0, t_e, t_p[i],  
                               t_s[i], t_lod_abs[i]);
    }
  }
}
