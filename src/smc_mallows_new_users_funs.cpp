#include <RcppArmadillo.h>
#include "missing_data.h"
#include "sample.h"
#include "parameterupdates.h"
#include "misc.h"
#include "partitionfuns.h"
#include "distances.h"
#include "classes.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


void reweight_new_users(
    SMCParameters& pars,
    const SMCAugmentation& aug,
    const SMCData& dat,
    const Rcpp::List& logz_list
){
  vec log_inc_wgt(pars.n_particles);
  for (int ii{}; ii < pars.n_particles; ++ii) {

    const double log_z_alpha = get_partition_function(
      dat.n_items, pars.alpha_samples(ii), logz_list, pars.metric
    );

    mat rankings = !aug.any_missing ? dat.rankings :
      aug.augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1),
        span(ii));

    double log_likelihood = -pars.alpha_samples(ii) / dat.n_items *
      rank_dist_sum(rankings, pars.rho_samples.col(ii), pars.metric,
                    dat.observation_frequency);

    log_inc_wgt(ii) = log_likelihood - dat.num_new_obs * log_z_alpha -
      log(aug.aug_prob(ii));
  }

  pars.norm_wgt = normalize_weights(log_inc_wgt);
  pars.effective_sample_size = (sum(pars.norm_wgt) * sum(pars.norm_wgt)) /
    sum(pars.norm_wgt % pars.norm_wgt);
}


