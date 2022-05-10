#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

cube initialize_rho(
  const int& N, const int& n_items, const int& d
){
  cube rho_samples(N, n_items, d, fill::zeros);
  for (uword i = 0; i < N; ++i) {
    uvec items_sample = randperm(n_items) + 1;
    for (uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }
  return rho_samples;
}

