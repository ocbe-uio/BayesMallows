#include <RcppArmadillo.h>
#include "distances.h"
#include "parameterupdates.h"
#include "sample.h"
#include "smc.h"
#include "partitionfuns.h"
#include "misc.h"
#include "smc_mallows_new_users.h"
#include "missing_data.h"

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
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const std::string aug_method = "uniform",
  const Rcpp::Nullable<arma::vec>& logz_estimate = R_NilValue,
  const Rcpp::Nullable<arma::vec>& cardinalities = R_NilValue,
  const std::string& metric = "footrule",
  const int& leap_size = 1,
  const Rcpp::Nullable<arma::cube> aug_init = R_NilValue,
  int num_obs = 0
) {

  int num_new_obs = new_rankings.n_cols;
  int n_users = rankings.n_cols;
  int n_items = rankings.n_rows;
  vec obs_freq = ones(n_users);
  vec aug_prob = ones(n_particles);
  bool any_missing = !is_finite(rankings);

  mat rho_samples(n_items, n_particles);
  vec alpha_samples = zeros(n_particles);
  double effective_sample_size;

  cube augmented_data;
  umat missing_indicator;


  if(any_missing){

    rankings.replace(datum::nan, 0);
    missing_indicator = conv_to<umat>::from(rankings);
    missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
    augmented_data = zeros(n_items, n_users, n_particles);

    for(int i{}; i < n_particles; i++) {
      augmented_data.slice(i) = rankings;
      initialize_missing_ranks(augmented_data.slice(i), missing_indicator);
    }

    if(aug_init.isNotNull()) {
      augmented_data(span::all, span(0, num_obs - 1), span::all) = Rcpp::as<cube>(aug_init);
    }
  } else {
    missing_indicator.reset();
  }

  num_obs += num_new_obs;
  mat new_observed_rankings, all_observed_rankings;
  if(!any_missing){
    new_observed_rankings = new_rankings;
    all_observed_rankings = rankings;
  }

  rho_samples = rho_init;
  alpha_samples = alpha_init;
  vec log_inc_wgt(n_particles, fill::zeros);

  if(any_missing){
    smc_mallows_new_users_augment_partial(
      augmented_data, aug_prob, rho_samples, alpha_samples, num_obs, num_new_obs,
      aug_method, missing_indicator, metric);
  }


  vec norm_wgt;
  smc_mallows_new_users_reweight(
    log_inc_wgt, effective_sample_size, norm_wgt, augmented_data, new_observed_rankings, rho_samples,
    alpha_samples, logz_estimate, cardinalities, num_obs, num_new_obs, aug_prob,
    any_missing, metric);

  smc_mallows_new_users_resample(
    rho_samples, alpha_samples, augmented_data, norm_wgt, num_obs,
    any_missing);

  for (int ii = 0; ii < n_particles; ++ii) {

    for (int kk = 0; kk < mcmc_steps; ++kk) {
      if(any_missing) {
        all_observed_rankings = augmented_data.slice(ii);
      }
      rho_samples.col(ii) = make_new_rho(rho_samples.col(ii), all_observed_rankings,
                      alpha_samples(ii), leap_size, metric, obs_freq);


      alpha_samples(ii) = update_alpha(alpha_samples(ii), all_observed_rankings,
                    obs_freq, rho_samples.col(ii), alpha_prop_sd, metric,
                    lambda, cardinalities, logz_estimate);

    }

    if(any_missing) {

      for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
        double log_hastings_correction = 0;

        uvec unranked_items = shuffle(find(missing_indicator.col(jj) == 1));
        if(aug_method == "uniform") {
          augmented_data(span::all, span(jj), span(ii)) = make_new_augmentation(
            augmented_data(span::all, span(jj), span(ii)),
            missing_indicator.col(jj),
            alpha_samples(ii),
            rho_samples.col(ii),
            metric
          );
        } else {

          Rcpp::List pprop = make_pseudo_proposal(
            unranked_items, augmented_data(span::all, span(jj), span(ii)),
            alpha_samples(ii), rho_samples.col(ii), metric, true
          );

          Rcpp::List bprop = make_pseudo_proposal(
            unranked_items, augmented_data(span::all, span(jj), span(ii)),
            alpha_samples(ii), rho_samples.col(ii), metric, false);
          double bprob = bprop["probability"];

          vec ar = pprop["proposal"];
          double prob = pprop["probability"];

          log_hastings_correction = -std::log(prob) + std::log(bprob);

          double ratio = -alpha_samples(ii) / n_items * (
            get_rank_distance(ar, rho_samples.col(ii), metric) -
              get_rank_distance(augmented_data(span::all, span(jj), span(ii)),
                                rho_samples.col(ii), metric)
          ) + log_hastings_correction;


          double u = std::log(randu<double>());

          if(ratio > u){
            augmented_data(span::all, span(jj), span(ii)) = ar;
          }
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
