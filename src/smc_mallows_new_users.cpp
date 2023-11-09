#include <RcppArmadillo.h>
#include "missing_data.h"
#include "smc_mallows_new_users.h"
#include "parameterupdates.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List  smc_mallows_new_users(
  arma::mat rankings,
  arma::mat new_rankings,
  const arma::mat rho_init,
  const arma::vec alpha_init,
  const int n_particles,
  const int mcmc_steps,
  const std::string aug_method,
  Rcpp::List logz_list,
  const std::string& metric,
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const int& leap_size = 1,
  const Rcpp::Nullable<arma::cube> aug_init = R_NilValue
) {

  int num_new_obs = new_rankings.n_cols;
  int n_users = rankings.n_cols;
  int n_items = rankings.n_rows;
  vec obs_freq = ones(n_users);
  vec aug_prob = ones(n_particles);
  bool any_missing = !is_finite(rankings);

  mat rho_samples = rho_init;
  vec alpha_samples = alpha_init;
  double effective_sample_size;

  cube augmented_data{};
  umat missing_indicator{};

  if(any_missing){
    set_up_missing(rankings, missing_indicator);
    augmented_data.set_size(n_items, n_users, n_particles);

    for(int i{}; i < n_particles; i++) {
      augmented_data.slice(i) = rankings;
      initialize_missing_ranks(augmented_data.slice(i), missing_indicator);
    }

    if(aug_init.isNotNull()) {
      augmented_data(span::all, span(0, rankings.n_cols - new_rankings.n_cols - 1), span::all) = Rcpp::as<cube>(aug_init);
    }
  }

  if(any_missing){
    smc_mallows_new_users_augment_partial(
      augmented_data, aug_prob, rho_samples, alpha_samples, num_new_obs,
      aug_method, missing_indicator, metric);
  }

  vec norm_wgt(n_particles);
  smc_mallows_new_users_reweight(
    effective_sample_size, norm_wgt, augmented_data, new_rankings, rho_samples,
    alpha_samples, logz_list, num_new_obs, aug_prob,
    any_missing, metric);

  smc_mallows_new_users_resample(
    rho_samples, alpha_samples, augmented_data, norm_wgt,
    any_missing);

  for (int ii = 0; ii < n_particles; ++ii) {

    for (int kk = 0; kk < mcmc_steps; ++kk) {
      if(any_missing) rankings = augmented_data.slice(ii);

      rho_samples.col(ii) =
        make_new_rho(rho_samples.col(ii), rankings,
                     alpha_samples(ii), leap_size, metric, obs_freq);

      alpha_samples(ii) = update_alpha(alpha_samples(ii), rankings,
                    obs_freq, rho_samples.col(ii), alpha_prop_sd, metric,
                    lambda, logz_list);

      if(any_missing) {
        int num_obs = rankings.n_cols;
        for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {

          augmented_data(span::all, span(jj), span(ii)) = make_new_augmentation(
            augmented_data(span::all, span(jj), span(ii)),
            missing_indicator.col(jj),
            alpha_samples(ii),
            rho_samples.col(ii),
            metric, aug_method == "pseudo"
          );
        }
      }

    }
  }

  // return the history of the particles and their values
  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = rho_samples,
    Rcpp::Named("alpha_samples") = alpha_samples,
    Rcpp::Named("augmented_rankings") = augmented_data,
    Rcpp::Named("ESS") = effective_sample_size
  );

  return particle_history;
}
