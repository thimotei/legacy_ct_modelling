vector onsets_lmpf(real inc_mean, real inc_sd, vector beta_im, vector beta_is,
                   matrix design, vector onset_avail, vector onset_time, 
                   vector onset_window, vector inf_time, array[] int ids) {
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
  }else{
    inc_sd_p = rep_vector(inc_sd, P);
  }

  // What is the probability of onset on observed day
  onset_upper = onset_time[ids] + inf_time[ids];
  onset_lower = fmax(rep_vector(0, nonsets), onset_upper - onset_window[ids]);
  for (i in 1:nonsets) {
    tar[ids[i]] = log_diff_exp(
      lognormal_lcdf(
        onset_upper[i] | inc_mean_p[ids[i]], inc_sd_p[ids[i]]
      ),
      lognormal_lcdf(
        onset_lower[i] | inc_mean_p[ids[i]], inc_sd_p[ids[i]]
      )
    );
  }

  return(tar);
}