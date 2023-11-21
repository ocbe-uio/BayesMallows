#include <RcppArmadillo.h>
#include "sample.h"
#include "partitionfuns.h"
#include "distances.h"
#include "misc.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


void update_dist_mat(mat& dist_mat, const mat& rankings,
                     const mat& rho_old, const std::string& metric,
                     const vec& observation_frequency){
  int n_clusters = dist_mat.n_cols;
  for(int i = 0; i < n_clusters; ++i)
    dist_mat.col(i) = rank_dist_vec(rankings, rho_old.col(i), metric, observation_frequency);
}






