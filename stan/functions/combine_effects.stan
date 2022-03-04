// Combine regression effects
// 
// Combines regression effects using a design matrix.
// 
// @param intercept The regression intercept
// 
// @param beta A Vector of effects. In general these should be specified
// on the unit scale as they may be rescaled (and hence pooled) using the
// beta_sd vector.
// 
// @param design The design matrix for the observations (rows) and effects
// columns.
// 
// @return A vector of linear predictions without error.
// 
// @author Sam Abbott
// @examples
// 
// # make example data
// intercept <- 1
// beta <- c(0.1, 0.2, 0.4)
// design <- t(matrix(c(1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1), 4, 4))
// 
// # Check effects  are combined as expected
// combine_effects(intercept, beta, design)
//
// Check function works with no effects
// combine_effects(intercept, as.double(c()), design)
vector combine_effects(real intercept, vector beta, matrix design) {
  int nobs = rows(design);
  int neffs = num_elements(beta);
  vector[neffs + 1] scaled_beta;
  if (neffs) {
    scaled_beta[1] = intercept;
    scaled_beta[2:(neffs+1)] = beta;
    return(design * scaled_beta);
  }else{
    return(rep_vector(intercept, nobs));
  }
}