#include <RcppArmadillo.h>
#include "missing_data.h"
#include "sample.h"
#include "parameterupdates.h"
#include "misc.h"
#include "partitionfuns.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]





void smc_mallows_new_users_resample(
    mat& rho_samples, vec& alpha_samples, cube& augmented_data,
    const vec& norm_wgt,
    const bool& partial
){
  int n_particles = rho_samples.n_cols;
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);
  rho_samples = rho_samples.cols(index);
  alpha_samples = alpha_samples.rows(index);

  if(partial) augmented_data = augmented_data.slices(index);
}

void smc_mallows_new_users_reweight(
    double& effective_sample_size,
    vec& norm_wgt,
    const cube& augmented_data,
    const mat& observed_rankings,
    const mat& rho_samples,
    const vec& alpha_samples,
    const Rcpp::List& logz_list,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& partial,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_cols;
  int n_items = rho_samples.n_rows;
  int num_obs = augmented_data.n_cols;
  vec log_inc_wgt(n_particles);
  for (int ii{}; ii < n_particles; ++ii) {

    const double log_z_alpha = get_partition_function(
      n_items, alpha_samples(ii), logz_list, metric
    );

    mat rankings = !partial ? observed_rankings : augmented_data(span::all,
      span(num_obs - num_new_obs, num_obs - 1),
      span(ii));

    double log_likelihood = -alpha_samples(ii) / n_items *
      rank_dist_sum(rankings, rho_samples.col(ii), metric, ones(rankings.n_cols));

    log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
  }

  norm_wgt = normalize_weights(log_inc_wgt);
  effective_sample_size = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);
}


