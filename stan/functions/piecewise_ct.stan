  vector piecewise_ct(vector t, real c0, real cp, real cs, 
                      real clod, real te, real tp, real ts,
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
    }else{
      lgrad = (clod - cp) / tlod;
      c_adj = cp;
    }

    for (i in 1:N) {
      if (adj_t[i] > 0) {
        if (adj_t[i] <= tp) {
          y[i] += adj_t[i] * ugrad;
        }else{
          adj_t[i] = adj_t[i] - tp;
          if (adj_t[i] <= ts) {
            y[i] = adj_t[i] * sgrad + cp;
          }else{
            adj_t[i] = adj_t[i] - ts;
            if (adj_t[i] <= tlod) {
              y[i] = adj_t[i] * lgrad + c_adj;
            }else{
              y[i] = clod;
            }
          }
        }
      }
    }
    return(y);
  }
  
  vector piecewise_ct_by_id(vector t, real c0, vector cp, vector cs, real clod, 
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
        t[t_start:t_end], c0, cp[k], cs[k], clod, te, tp[k], ts[k], tlod[k],
        switch
      );
    }

    return(ct);
  }
