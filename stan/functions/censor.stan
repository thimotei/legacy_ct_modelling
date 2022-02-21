real censor(real mu, real bound) {
  real y = mu > bound ? bound : mu;
  return y;
}
