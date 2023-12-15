#include <RcppArmadillo.h>
#include "classes.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List  run_smc(
  Rcpp::List data,
  Rcpp::List new_data,
  Rcpp::List model_options,
  Rcpp::List smc_options,
  Rcpp::List compute_options,
  Rcpp::List priors,
  Rcpp::List initial_values,
  Rcpp::Nullable<arma::mat> pfun_values,
  Rcpp::Nullable<arma::mat> pfun_estimate) {

  SMCData dat{data, new_data};
  SMCParameters pars{model_options, smc_options, compute_options, initial_values};
  Priors pris{priors};
  SMCAugmentation aug{dat, smc_options, initial_values, pars.n_particles};
  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  aug.reweight(pars, dat, pfun);
  uvec index = pars.draw_resampling_index();
  pars.resample(index);
  aug.resample(index);

  for (size_t ii{}; ii < pars.n_particles; ++ii) {
    for (size_t kk{}; kk < pars.mcmc_steps; ++kk) {
      aug.update_data(ii, dat);
      pars.update_rho(ii, dat);
      pars.update_alpha(ii, dat, pfun, pris);
      aug.update_missing_ranks(ii, dat, pars);
    }
  }

  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = pars.rho_samples,
    Rcpp::Named("alpha_samples") = pars.alpha_samples,
    Rcpp::Named("augmented_rankings") = aug.augmented_data
  );

  return particle_history;
}
