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





vec update_wcd(const uvec& current_cluster_assignment,
                     const mat& dist_mat){
  int n_clusters = dist_mat.n_cols;
  vec wcd(n_clusters);

  uvec inds = regspace<uvec>(0, n_clusters - 1);
  for(int i = 0; i < n_clusters; ++i){
    mat dist_vec = dist_mat.submat(find(current_cluster_assignment == i), inds.subvec(i, i));
    wcd(i) = accu(dist_vec);
  }

  return(wcd);
}
