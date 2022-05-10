#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void smc_mallows_new_users_reweight(
    vec& log_inc_wgt,
    rowvec& ESS_vec,
    vec& norm_wgt,
    const cube& aug_rankings,
    const cube& rho_samples,
    const double& alpha,
    const mat& alpha_samples,
    const int& N,
    const int& tt,
    const int& n_items,
    const Rcpp::Nullable<vec> logz_estimate,
    const std::string& metric,
    const int& num_obs,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& augment_alpha
){
  for (int ii{}; ii < N; ++ii) {
    Rcpp::Nullable<vec> cardinalities = R_NilValue;

    rowvec rho_samples_ii = \
      rho_samples(span(ii), span::all, span(tt + 1));

    double alpha_used = augment_alpha ? alpha_samples(ii, tt + 1) : alpha;

    /* Calculating log_z_alpha and log_likelihood ----------- */
    const double log_z_alpha = get_partition_function(
      n_items, alpha_used, cardinalities, logz_estimate, metric
    );

    mat new_observed_rankings;
    new_observed_rankings = aug_rankings(span(num_obs - num_new_obs, num_obs - 1), span::all, span(ii));
    double log_likelihood = get_exponent_sum(                          \
      alpha_used, rho_samples_ii.t(), n_items, new_observed_rankings, metric\
    );
    log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
  }

  /* normalise weights ------------------------------------ */
  const double maxw = max(log_inc_wgt);
  const vec w = exp(log_inc_wgt - maxw);
  norm_wgt = w / sum(w);

  ESS_vec(tt) = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);

}

