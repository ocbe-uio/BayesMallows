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
  if(dat.num_new_obs == 0) return;
  vec log_inc_wgt(pars.n_particles);
  for (size_t particle{}; particle < pars.n_particles; ++particle) {

    const double log_z_alpha = get_partition_function(
      dat.n_items, pars.alpha_samples(particle), logz_list, pars.metric
    );

    mat new_rankings = !aug.any_missing ? dat.new_rankings :
      aug.augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1),
        span(particle));

    double log_likelihood = -pars.alpha_samples(particle) / dat.n_items *
      rank_dist_sum(new_rankings, pars.rho_samples.col(particle), pars.metric,
                    dat.observation_frequency(span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1)));

    log_inc_wgt(particle) = log_likelihood - dat.num_new_obs * log_z_alpha -
      log(aug.aug_prob(particle));
  }

  pars.norm_wgt = normalize_weights(log_inc_wgt);
  pars.effective_sample_size = 1 / std::pow(norm(pars.norm_wgt, 2), 2);
}


