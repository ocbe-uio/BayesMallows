#include <RcppArmadillo.h>
#include "missing_data.h"
#include "smc_mallows_new_users.h"
#include "parameterupdates.h"
#include "classes.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List  run_smc(
  Rcpp::List data,
  Rcpp::List new_data,
  Rcpp::List smc_options,
  Rcpp::List compute_options,
  Rcpp::List initial_values,
  Rcpp::List logz_list,
  const std::string& metric,
  const double lambda = 0.1
) {

  SMCData dat{data, new_data, Rcpp::List{}};
  SMCOptions smc_opt{smc_options};
  SMCParameters pars{compute_options, initial_values};

  vec aug_prob = ones(smc_opt.n_particles);
  bool any_missing = !is_finite(dat.rankings);

  double effective_sample_size;

  cube augmented_data{};
  umat missing_indicator{};

  if(any_missing){
    set_up_missing(dat.rankings, missing_indicator);
    augmented_data.set_size(dat.n_items, dat.n_assessors, smc_opt.n_particles);

    for(int i{}; i < smc_opt.n_particles; i++) {
      augmented_data.slice(i) = dat.rankings;
      initialize_missing_ranks(augmented_data.slice(i), missing_indicator);
    }

    Rcpp::Nullable<arma::cube> aug_init = initial_values["aug_init"];
    if(aug_init.isNotNull()) {
      augmented_data(span::all, span(0, dat.rankings.n_cols - dat.new_rankings.n_cols - 1), span::all) = Rcpp::as<cube>(aug_init);
    }
  }

  if(any_missing){
    smc_mallows_new_users_augment_partial(
      augmented_data, aug_prob, pars.rho_samples, pars.alpha_samples, dat.num_new_obs,
      smc_opt.aug_method, missing_indicator, metric);
  }

  vec norm_wgt(smc_opt.n_particles);
  smc_mallows_new_users_reweight(
    effective_sample_size, norm_wgt, augmented_data, dat.new_rankings, pars.rho_samples,
    pars.alpha_samples, logz_list, dat.num_new_obs, aug_prob,
    any_missing, metric);

  smc_mallows_new_users_resample(
    pars.rho_samples, pars.alpha_samples, augmented_data, norm_wgt,
    any_missing);

  for (int ii = 0; ii < smc_opt.n_particles; ++ii) {

    for (int kk = 0; kk < smc_opt.mcmc_steps; ++kk) {
      if(any_missing) dat.rankings = augmented_data.slice(ii);

      pars.rho_samples.col(ii) =
        make_new_rho(pars.rho_samples.col(ii), dat.rankings,
                     pars.alpha_samples(ii), pars.leap_size, metric, dat.observation_frequency);

      pars.alpha_samples(ii) = update_alpha(pars.alpha_samples(ii), dat.rankings,
                    dat.observation_frequency, pars.rho_samples.col(ii), pars.alpha_prop_sd, metric,
                    lambda, logz_list);

      if(any_missing) {
        int num_obs = dat.rankings.n_cols;
        for (int jj = num_obs - dat.num_new_obs; jj < num_obs; ++jj) {

          augmented_data(span::all, span(jj), span(ii)) = make_new_augmentation(
            augmented_data(span::all, span(jj), span(ii)),
            missing_indicator.col(jj),
            pars.alpha_samples(ii),
            pars.rho_samples.col(ii),
            metric, smc_opt.aug_method == "pseudo"
          );
        }
      }

    }
  }

  // return the history of the particles and their values
  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = pars.rho_samples,
    Rcpp::Named("alpha_samples") = pars.alpha_samples,
    Rcpp::Named("augmented_rankings") = augmented_data,
    Rcpp::Named("ESS") = effective_sample_size
  );

  return particle_history;
}
