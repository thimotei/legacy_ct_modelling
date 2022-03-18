vector onsets_lmpf(real inc_mean, real inc_sd, vector beta_im, vector beta_is,
                   matrix design, vector onset_avail, vector onset_time, 
                   vector inf_time) {
  int P = num_elements(onset_avail);
  vector[P] inc_mean_p;
  vector[P] inc_sd_p;
  vector[P] tar = rep_vector(0, P);
  int adj_inc_sd = num_elements(beta_is);

  inc_mean_p = combine_effects(inc_mean, beta_im, design);
  if (adj_inc_sd) {
    inc_sd_p = exp(combine_effects(log(inc_sd), beta_is, design));
  }else{
    inc_sd_p = rep_vector(inc_sd, P);
  }

   // Component of likelihood for time of exposure
   for(j in 1:P) {
   // likelihood for time of exposure using the CDF of the incubation period
   // and known symptom onset day
   // What is the probability onsets on observed day
    if (onset_avail[j]) {
      real onset_from_inf = onset_time[j] + inf_time[j];
      real onset_window = max({0, onset_from_inf - 1});
        tar[j] = log_diff_exp(
          lognormal_lcdf(onset_from_inf | inc_mean_p[j], inc_sd_p[j]),
          lognormal_lcdf(onset_window | inc_mean_p[j], inc_sd_p[j])
        );
    }
  }
  return(tar);
}