#include <RcppArmadillo.h>
#include "sample.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void smc_mallows_new_users_resample(
    cube& rho_samples,
    mat& alpha_samples,
    cube& aug_rankings,
    const int& N,
    const vec& norm_wgt,
    const int& tt,
    const int& num_obs,
    const bool& augment_alpha,
    const bool& partial
){
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, N - 1), N, true, norm_wgt);
  uvec tt_vec;
  tt_vec = tt;

  /* Replacing tt + 1 slice on rho_samples ---------------- */
  mat rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
  rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
  rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

  /* Replacing tt + 1 column on alpha_samples ------------- */
  if(augment_alpha){
    alpha_samples.col(tt + 1) =  alpha_samples.submat(index, tt_vec + 1);
  }

  if(partial){
    cube aug_rankings_index = aug_rankings.slices(index);
    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(span(0, num_obs - 1), span::all, span::all);
  }

}

