// #include functions/piecewise_ct.stan
// #include functions/combine_effects.stan
// #include functions/onsets_lmpf.stan

functions{
  vector combine_effects(real intercept, vector beta, matrix design) {
    
    int nobs = rows(design);
    int neffs = num_elements(beta);
    vector[neffs + 1] scaled_beta; // +1 dimensions to include the single intercept
    
    if (neffs) {
      
      scaled_beta[1] = intercept;
      scaled_beta[2:(neffs+1)] = beta;
      
      return(design * scaled_beta);
      
    } else {
      
      return(rep_vector(intercept, nobs));
      
    }
  }
  
  vector onsets_lmpf(
    real inc_mean, real inc_sd, vector beta_im, vector beta_is, matrix design,
    vector onset_avail, vector onset_time, vector onset_window, vector inf_time,
    array[] int ids) {
      
      int P = num_elements(onset_avail);
      vector[P] inc_mean_p;
      vector[P] inc_sd_p;
      vector[P] tar = rep_vector(0, P);
      int adj_inc_sd = num_elements(beta_is);
      int nonsets = num_elements(ids);
      vector[nonsets] onset_upper;
      vector[nonsets] onset_lower;
      
      inc_mean_p = combine_effects(inc_mean, beta_im, design);
      
      if (adj_inc_sd) {
        inc_sd_p = exp(combine_effects(log(inc_sd), beta_is, design));
      } else {
        inc_sd_p = rep_vector(inc_sd, P);
      }
      
      // What is the probability of onset on observed day
      onset_upper = onset_time[ids] + inf_time[ids];
      onset_lower = fmax(rep_vector(0, nonsets), 
                         onset_upper - onset_window[ids]);
      
      for (i in 1:nonsets) {
        
        tar[ids[i]] = log_diff_exp(
          lognormal_lcdf(onset_upper[i] | inc_mean_p[ids[i]], inc_sd_p[ids[i]]),
          lognormal_lcdf(onset_lower[i] | inc_mean_p[ids[i]], inc_sd_p[ids[i]])
          );
          
      }
      
      return(tar);
    }
    
    vector piecewise_ct(
      vector t, real c0, real cp, real cs, real clod, real te, real tp, real ts,
      real tlod, int switch) {
        
        int N = num_elements(t);
        vector[N] y = rep_vector(c0, N);
        vector[N] adj_t = t - te;
        real ugrad = (cp - c0) / tp;
        real sgrad;
        real lgrad;
        real c_adj;
        
        if (switch) {
          
          sgrad = (cs - cp) / ts;
          lgrad = (clod - cs) / tlod;
          c_adj = cs;
          
        } else {
          
          lgrad = (clod - cp) / tlod;
          c_adj = cp;
          
        }
        
        for (i in 1:N) {
          if (adj_t[i] > 0) {
            if (adj_t[i] <= tp) {
              
              y[i] += adj_t[i] * ugrad;
              
            } else {
              
              adj_t[i] = adj_t[i] - tp;
              
              if (adj_t[i] <= ts) {
                
                y[i] = adj_t[i] * sgrad + cp;
                
              } else {
                
                adj_t[i] = adj_t[i] - ts;
                if (adj_t[i] <= tlod) {
                  
                  y[i] = adj_t[i] * lgrad + c_adj;
                  
                } else {
                  
                  y[i] = clod;
                  
                }
              }
            }
          }
        }
        
        return(y);
      }
    
    vector piecewise_ct_by_id(
      vector t, real c0, vector cp, vector cs, real clod, 
      real te, vector tp, vector ts, vector tlod,
      array[] int id, array[] int tests_per_id,
      array[] int cum_tests_per_id, int switch) {
        int P = num_elements(cp);
        int N = num_elements(t);
        vector[N] ct;
        
        for(k in 1:P) {
          int t_start = cum_tests_per_id[k] - tests_per_id[k] + 1;
          int t_end = cum_tests_per_id[k];
          
          ct[t_start:t_end] = piecewise_ct(
            t[t_start:t_end], c0, cp[k], cs[k], clod, 
            te, tp[k], ts[k], tlod[k], switch);
        }
        
        return(ct);
      }
      
      
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
  int ncensored; // Number of censored tests
  array[ncensored] int censored; // Which tests have been censored
  int nuncensored; // Number of uncensored tests
  array[nuncensored] int uncensored; // Which tests have not been censored
  array[N] int uncensored_by_test; // Uncensored by test
  vector[N] day_rel; // day of test (integer)
  vector<lower = 0>[N] ct_value; // Ct value of test
  int any_onsets; // Are there any symptom onsets
  int nonsets; // Number of onsets
  vector[P] onset_avail; // Onsets available per ID
  vector[P] onset_time; // Time of onset per ID
  vector[P] onset_window; // Window in which onsets could have occurred per ID
  array[nonsets] int ids_with_onsets; // IDs that have onsets
  int K; //Number of parameters with individual level variation
  int switch; //Should a secondary breakpoint in the CT curve be modelled
  int ind_var_m; // Should inividual variation be modelled
  real ind_var_sd; // Standard deviation of inividual variation to be modelled
  int ind_corr; // Should individual variation be modelled with correlation
  real lkj_prior; // LKJ prior for individual level variation
  int preds; // Number of predictors
  real preds_sd; // Standard deviation of predictor coeffs
  matrix[P, preds + 1] design; //Design matrix
  int adj_t_p; // Should time at peak be adjusted
  int adj_t_s; // Should time at switch be adjusted
  int adj_t_lod; // Should time at LOD be adjusted
  int adj_c_p; // Should Ct at peak be adjusted
  int adj_c_s; // Should Ct at switch be adjusted
  int adj_inc_mean;// Should incubation period mean be adjusted
  int adj_inc_sd; // Should incubation period standard deviation be adjusted
  int adj_ct; // Should Cts be adjusted
  int ct_preds; // Number of predictors for CT adjustment
  real ct_preds_sd; // Standard deviation of CT predictor coeffs
  matrix[N, ct_preds + 1] ct_design; // Design matrix for CT adjustment
  int likelihood;
  int output_loglik;
}

transformed data {
  
  vector[P] t_inf_bound;
  vector[61] sim_times;
  
  for (i in 1:P) {
    t_inf_bound[i] = max({-onset_time[i], 0});
  }
  for (i in 0:60) {
    sim_times[i + 1] = i;
  }
  
}

parameters {
  
  vector<lower = t_inf_bound>[P] t_inf; // Inferred time of infection
  array[any_onsets] real inc_mean; //Incubation period mean
  array[any_onsets] real<lower = 0> inc_sd; //Incubation period sd
  real<lower = c_lod> c_0;   // Ct value before detection
  
  // Cholesky factored correlation matrix
  cholesky_factor_corr[ind_corr ? K : 0] L_Omega;
  real c_p_mean; // Ct value at peak
  array[switch] real c_s_mean; // Ct value at switch
  real t_p_mean; // Timing of peak
  array[switch] real t_s_mean; // Timing of switch
  real t_lod_mean; // Time viral load hits lower limit of detection
  vector<lower = 0>[ind_var_m ? K : 0] ind_var; // SD of individual variation
  matrix[ind_var_m ? K : 0, P] ind_eta; // Individual level variation
  real<lower = 0> sigma; // Variance parameter for observation model
  
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
  
  // individual-level Ct dynamics parameters
  vector[P] t_p, t_s, t_lod, c_p, c_s, t_lod_abs;
  vector[N] inf_rel, exp_ct, adj_exp_ct;
  
  // likelihoods for onset and incubation periods
  array[nonsets] real onsets_star;
  vector[nonsets ? P :0] onsets_log_lik; 
  
  {
    matrix[P, K] eta;
    
    if (ind_corr) {
      // Cholesky factor of the covariance matrix
      matrix[K, K] L;
      L = diag_pre_multiply(ind_var, L_Omega);
      
      // Calculate per-infection correlated effects
      eta = (L * ind_eta)';
      
    } else {
      if (ind_var_m) {
        
        // All effects are independent
        for (i in 1:K) {
          
          eta[1:P, i] = to_vector(ind_eta[i, 1:P]) * ind_var[i];
          
        }
      } else {
        
        // No infection level differences
        eta = rep_matrix(0, P, K);
      }
      
    }
    
    // Combine effects for each CT parameter and transform to required scale
    t_p = exp(combine_effects(t_p_mean, beta_t_p, design) + eta[, 1]);
    t_lod = exp(combine_effects(t_lod_mean, beta_t_lod, design) + eta[, 2]);
    c_p = inv_logit(combine_effects(c_p_mean, beta_c_p, design) + eta[, 3]);
    
    // Optional effects if a second breakpoint is used
    if (switch) {
      
      t_s = exp(combine_effects(t_s_mean[1], beta_t_s, design) + eta[, 4]);
      c_s = c_0 * inv_logit(
        combine_effects(c_s_mean[1], beta_c_s, design) + eta[, 5]
        );
        c_p = c_s .* c_p;
        
    } else {
      
      c_s = rep_vector(0.0, P);
      t_s = rep_vector(0.0, P);
      c_p = c_0 * c_p;
      
    }
  }
  // Make times absolute
  t_lod_abs = t_p + t_s + t_lod;
  
  // Adjust observed times since first test to be time since infection
  inf_rel = day_rel + t_inf[id];

  // Expected ct value given viral load parameters
  exp_ct = piecewise_ct_by_id(
    inf_rel, c_0, c_p, c_s, c_0, t_e, t_p, t_s, t_lod, id,
    tests_per_id, cum_tests_per_id, switch
  );

  // Shift and scale ct values
  adj_exp_ct = combine_effects(0, beta_ct_shift, ct_design) +
    exp(combine_effects(0, beta_ct_scale, ct_design)) .* exp_ct;

  // Model symptom onset likelihood: see onsets_lmpf.stan
  if (any_onsets) {
    
    vector[P] onsets_ttar;

    onsets_ttar = onsets_lmpf(
      inc_mean[1], inc_sd[1], beta_inc_mean, beta_inc_sd, design, onset_avail,
      onset_time, onset_window, t_inf, ids_with_onsets
    );
    
    onsets_star[1] = sum(onsets_ttar);
    if (output_loglik) {
      onsets_log_lik = onsets_ttar;
    }
  }
}

model {
  // Prior over possible infection times relative to first
  // positive test or symtom onset.
  // Assumes that the first positive test is not a false positive.
  t_inf ~ normal(t_inf_bound + 5, 5); 
  target += -normal_lccdf(t_inf_bound | t_inf, 5);

  // CT piecewise linear intercept parameters
  c_0 ~ normal(c_lod + 10, 2) T[c_lod, ];
  c_p_mean ~ normal(0, 1); //mean at 50% of switch value
  t_p_mean ~ normal(1.6, 0.5); //mean at log(5)
  t_lod_mean ~ normal(2.3, 0.5); //mean at log(10) + peak + scale timing 
  if (switch) {
    c_s_mean ~ normal(0, 1); //mean at 50% of maximum ct
    t_s_mean ~ normal(1.61, 0.5); //mean at log(5) + peak timing
  }

  // Individual level variation
  if (ind_var_m) {
    to_vector(ind_eta) ~ std_normal();
    ind_var ~ normal(0, ind_var_sd);
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

  if (any_onsets) {
    // Priors on the incubation period
    inc_mean ~ normal(lmean[1], lmean[2]);
    inc_sd[1] ~ normal(lsd[1], lsd[2]) T[0, ];
 
    if (likelihood) {
      // Component of likelihood for symptom onsets see onsets_lpmf.stan
      target += onsets_star[1];
    }
  }

  if (likelihood) {
    // Component of likelihood for expected ct values
    // If non-censored: P(observed ct | expected ct)
    ct_value[uncensored] ~ normal(adj_exp_ct[uncensored], sigma);
    // If censored: P(expected ct >= censored ct)
    target += normal_lccdf(c_lod | adj_exp_ct[censored], sigma);
    // All CTs are truncated above 0. P(0 <= expected ct)
    target += -normal_lccdf(0 | adj_exp_ct, sigma);
  }
}

generated quantities {
  vector[N] sim_ct;
  matrix[ind_corr ? K : 0, ind_corr ? K : 0] correlation;
  vector[output_loglik ? N : 0] log_lik;
  if (ind_corr) {
    correlation = L_Omega * L_Omega';
  }
  // Posterior predictions
  sim_ct = to_vector(normal_rng(adj_exp_ct, sigma));
  sim_ct = fmin(sim_ct, c_lod);
  // Output by infection log-likelihood
  if (output_loglik) {
    log_lik = rep_vector(0, N);
    for(i in 1:N){
      if (uncensored_by_test[i] == 1) {
        log_lik[i] = normal_lpdf(ct_value[i] | adj_exp_ct[i], sigma) - normal_lccdf(0 | adj_exp_ct[i], sigma);
      } else {
        log_lik[i] = normal_lccdf(c_lod | adj_exp_ct[i], sigma) - normal_lccdf(0 | adj_exp_ct[i], sigma);
      }
    }
  }
}