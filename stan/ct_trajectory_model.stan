functions{ 
#include functions/piecewise_ct.stan
#include functions/combine_effects.stan
#include functions/truncated_normal_rng.stan
#include functions/censor.stan
}

data {
  int P; // number of patients
  int N; // number of tests
  array[N] int id; // id of person
  array[P] int tests_per_id; // Tests per ID
  array[P] int cum_tests_per_id; // Cumulative tests per id
  real c_lod; // Ct value at limit of detection 
  real t_e; 
  array[2] real lmean; // mean of incubation period used (+ sd)
  array[2] real lsd; // standard deviation of incubation period used (+ sd)
  array[N] int pcr_res; // boolean test result
  vector[N] day_rel; // day of test (integer)
  vector[N] ct_value; // Ct value of test
  int any_onsets;
  vector[P] onset_avail;
  vector[P] onset_time;
  int K; //Number of parameters with individual level variation
  int switch; //Should a secondary breakpoint in the CT curve be modelled
  int ind_var_m; // Should inividual variation be modelled
  int ind_corr; // Should individual variation be modelled with correlation
  real lkj_prior; // LKJ prior for individual level variation
  int preds; // Number of predictors
  real preds_sd; // Standard deviation of predictor coeffs
  matrix[P, preds + 1] design; //Design matrix
  int adj_t_p; // Should time at peak be adjusted
  int adj_t_s; // Should time at switch be adjusted
  int adj_t_lod; // Should time at LOD be adjusted
  int adj_c_p; // Should CT at peak be adjusted
  int adj_c_s; // Should CT at switch be adjusted
  int adj_inc_mean; // Should incubation period mean be adjusted
  int adj_inc_sd; // Should incubation period standard deviation be adjusted
  int adj_ct; // Should cts be adjusted
  int ct_preds; // Number of predictors for CT adjustment
  real ct_preds_sd; // Standard deviation of CT predictor coeffs
  matrix[N, ct_preds + 1] ct_design; // Design matrix for CT adjustment
  int likelihood;
  int output_loglik;
}

transformed data {
  vector[P] T_e_bound;
  vector[61] sim_times;
  for (i in 1:P) {
    T_e_bound[i] = max({-onset_time[i], 0});
  }
  for (i in 0:60) {
    sim_times[i + 1] = i;
  }
}

parameters {
  vector<lower = T_e_bound>[P] T_e; // Inferred time of infection
  array[any_onsets] real inc_mean; //Incubation period mean
  array[any_onsets] real<lower = 0> inc_sd; //Incubation period sd
  real<lower = c_lod> c_0;   // Ct value before detection
  // Cholesky_factored correlation matrix
  cholesky_factor_corr[ind_corr ? K : 0] L_Omega;
  real c_p_mean; // Ct value of viral load p
  array[switch] real c_s_mean; // Ct value at s
  real t_p_mean; // Timing of peak
  array[switch] real t_s_mean; // Timing of switch
  real t_lod_mean; // Time viral load hits lower limit of detection
  vector<lower = 0>[ind_var_m ? K : 0] ind_var; // SD of individual variation
  matrix[ind_var_m ? K : 0, P] ind_eta; // Individual level variation
  real<lower = 0> sigma; // Variance parameter for oobservation model
  // Coefficients
  vector[preds && adj_t_p ? preds : 0] beta_t_p;
  vector[preds && adj_t_s ? preds : 0] beta_t_s;
  vector[preds && adj_t_lod ? preds : 0] beta_t_lod;
  vector[preds && adj_c_s ? preds : 0] beta_c_s;
  vector[preds && adj_c_p ? preds : 0] beta_c_p;
  vector[preds && adj_inc_mean ? preds : 0] beta_inc_mean;
  vector[preds && adj_inc_sd ? preds : 0] beta_inc_sd;
  vector[ct_preds && adj_ct ? ct_preds : 0] beta_ct_shift;
  vector[ct_preds && adj_ct ? ct_preds : 0] beta_ct_scale;
}

transformed parameters {
  vector[P] t_p; vector[P] t_s; vector[P] t_lod;
  vector[P] c_p; vector[P] c_s;
  vector[P] t_lod_abs; vector[N] t_inf;
  vector[N] exp_ct; vector[N] adj_exp_ct;
{
  matrix[P, K] eta;
  if (ind_corr) {
    // Cholesky factor of the covariance matrix
    matrix[K, K] L;
    L = diag_pre_multiply(ind_var, L_Omega);

    // Calculate per-person correlated effects
    eta = (L * ind_eta)';
  }else{
    if (ind_var_m) {
      for (i in 1:K) {
        eta[1:P, i] = to_vector(ind_eta[i, 1:P]) * ind_var[i];
      }
    }else{
      eta = rep_matrix(0, P, K);
    }
  }

  // Combine effects for each CT parameter and transform to required scale
  t_p = exp(combine_effects(t_p_mean, beta_t_p, design)+ eta[, 1]);
  t_lod = exp(combine_effects(t_lod_mean, beta_t_lod, design) + eta[, 2]);
  c_p = inv_logit(combine_effects(c_p_mean, beta_c_p, design)+ eta[, 3]);
  // Optional effects if a second breakpoint is used
  if (switch) {
    t_s = exp(combine_effects(t_s_mean[1], beta_t_s, design) + eta[, 4]);
    c_s = c_0 * inv_logit(
      combine_effects(c_s_mean[1], beta_c_s, design) + eta[, 5]
    );
    c_p = c_s .* c_p;
  }else{
    c_s = rep_vector(0.0, P);
    t_s = rep_vector(0.0, P);
    c_p = c_0 * c_p;
  }
}
  // Make times absolute
  t_lod_abs = t_p + t_s + t_lod;
  // Adjust observed times since first test to be time since infection
  t_inf = day_rel + T_e[id];

  // Expected ct value given viral load parameters
  exp_ct = piecewise_ct_by_id(
    t_inf, c_0, c_p, c_s, c_0, t_e, t_p, t_s, t_lod_abs, id,
    tests_per_id, cum_tests_per_id, switch
  );

  // Shift and scale ct values
  adj_exp_ct = combine_effects(0, beta_ct_shift, ct_design) +
    exp(combine_effects(0, beta_ct_scale, ct_design)) .* exp_ct;
}

model {
  // Prior over possible infection times relative to first
  // positive test or symtom onset.
  // Assumes that the first positive test is not a false positive.
  for (i in 1:P) {
    T_e[i] ~ normal(T_e_bound[i] + 5, 5) T[T_e_bound[i],];
  }

  // CT piecewise linear intercept parameters
  c_0 ~ normal(c_lod + 10, 5) T[0, ];
  c_p_mean ~ normal(0, 1); //mean at 50% of switch value
  t_p_mean ~ normal(1.61, 0.5); //mean at log(5)
  t_lod_mean ~ normal(2.3, 0.5); //mean at log(10) + peak + scale timing 
  if (switch) {
    c_s_mean ~ normal(0, 1); //mean at 50% of maximum ct
    t_s_mean ~ normal(1.61, 0.5); //mean at log(5) + peak timing
  }

  // Individual level variation
  if (ind_var_m) {
    to_vector(ind_eta) ~ std_normal();
    ind_var ~ normal(0, 0.25);
  }

  // LKJ prior on correlation between individual level dynamics
  if (ind_corr) {
    L_Omega ~ lkj_corr_cholesky(lkj_prior);
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
      beta_t_lod ~ normal(0, preds_sd);
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
    if (adj_inc_sd) {
      beta_inc_sd ~ normal(0, preds_sd);
    } 
  }
  if (ct_preds && adj_ct) {
    beta_ct_shift ~ normal(0, ct_preds_sd);
    beta_ct_scale ~ normal(0, ct_preds_sd);
  }

  if (any_onsets && likelihood) {
    vector[P] inc_mean_p;
    vector[P] inc_sd_p;
    // Priors on the incubation period
    inc_mean[1] ~ normal(lmean[1], lmean[2]);
    inc_mean_p = combine_effects(inc_mean[1], beta_inc_mean, design);
    inc_sd[1] ~ normal(lsd[1], lsd[2]) T[0, ];
    if (adj_inc_sd) {
      inc_sd_p = exp(combine_effects(log(inc_sd[1]), beta_inc_sd, design));
    }else{
      inc_sd_p = rep_vector(inc_sd[1], P);
    }

   // Component of likelihood for time of exposure
   for(j in 1:P) {
   // likelihood for time of exposure using the CDF of the incubation period
   // and known symptom onset day
   // What is the probability onsets on observed day
    if (onset_avail[j]) {
      real onset_from_inf = onset_time[j] + T_e[j];
      real onset_window = max({0, onset_from_inf - 1});
        target += log_diff_exp(
          lognormal_lcdf(onset_from_inf | inc_mean_p[j], inc_sd_p[j]),
          lognormal_lcdf(onset_window | inc_mean_p[j], inc_sd_p[j])
        );
     }
    }
  }

  if (likelihood) {
    // Component of likelihood for expected ct values
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
  vector[N] sim_ct;
  matrix[ind_corr ? K : 0, ind_corr ? K : 0] correlation;
  if (ind_corr) {
    correlation = L_Omega * L_Omega';
  }
  for (i in 1:N) {
    sim_ct[i] = truncated_normal_rng(adj_exp_ct[i], sigma, 0, c_0);
    sim_ct[i] = censor(sim_ct[i], c_lod);
  }
  if (output_loglik) {

  }
}
