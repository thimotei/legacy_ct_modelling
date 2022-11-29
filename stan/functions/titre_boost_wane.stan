vector titre_boost_wane(vector t, real tb, real tp, real ts, 
                        real m1, real m2, real m3) {

  int N = num_elements(t);
  vector[N] y = rep_vector(0, N);
  
  for (i in 1:N) {
    if (t[i] > 0 && t[i] <= tb) {
      y[i] = 0;
    }
    else if(t > tb && t <= tp + tb) {
      y[i] = m1 * (t - tb);
    }
    else if(t > tp + tb && t <= tp + tb + ts) {
      y[i] = m2*(t - tb - tp) + m1*tp;
    }
    else if(t > tp + tb + ts) {
      y[i] = m3*(t - tb - tp - ts) + m1*tp + m2*ts;
    }
  }
    
  return(y);
}
