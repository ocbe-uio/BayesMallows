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
  Rcpp::List model_options,
  Rcpp::List smc_options,
  Rcpp::List compute_options,
  Rcpp::List priors,
  Rcpp::List initial_values,
  Rcpp::List logz_list
) {

  SMCData dat{data, new_data};
  SMCParameters pars{model_options, smc_options, compute_options, initial_values};
  Priors pris{priors};
  SMCAugmentation aug{dat, smc_options, initial_values, pars.n_particles};

  aug.augment_partial(pars, dat);
  reweight_new_users(pars, aug, dat, logz_list);

  uvec index = pars.draw_resampling_index();
  pars.resample(index);
  aug.resample(index);


  for (int ii = 0; ii < pars.n_particles; ++ii) {

    for (int kk = 0; kk < pars.mcmc_steps; ++kk) {
      if(aug.any_missing) dat.rankings = aug.augmented_data.slice(ii);

      pars.rho_samples.col(ii) =
        make_new_rho(pars.rho_samples.col(ii), dat.rankings,
                     pars.alpha_samples(ii), pars.leap_size, pars.metric, dat.observation_frequency);

      pars.alpha_samples(ii) = update_alpha(pars.alpha_samples(ii), dat.rankings,
                    dat.observation_frequency, pars.rho_samples.col(ii), pars.alpha_prop_sd, pars.metric,
                    pris.lambda, logz_list);

      if(aug.any_missing) {
        int num_obs = dat.rankings.n_cols;
        for (int jj = num_obs - dat.num_new_obs; jj < num_obs; ++jj) {

          aug.augmented_data(span::all, span(jj), span(ii)) = make_new_augmentation(
            aug.augmented_data(span::all, span(jj), span(ii)),
            aug.missing_indicator.col(jj),
            pars.alpha_samples(ii),
            pars.rho_samples.col(ii),
            pars.metric, aug.aug_method == "pseudo"
          );
        }
      }

    }
  }

  // return the history of the particles and their values
  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = pars.rho_samples,
    Rcpp::Named("alpha_samples") = pars.alpha_samples,
    Rcpp::Named("augmented_rankings") = aug.augmented_data,
    Rcpp::Named("ESS") = pars.effective_sample_size
  );

  return particle_history;
}
