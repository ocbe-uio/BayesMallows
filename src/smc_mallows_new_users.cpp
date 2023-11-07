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
  arma::mat rho_init,
  arma::vec alpha_init,
  const std::string& type,
  const int& n_particles,
  const int& mcmc_steps,
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const std::string& aug_method = "uniform",
  const Rcpp::Nullable<arma::vec>& logz_estimate = R_NilValue,
  const Rcpp::Nullable<arma::vec>& cardinalities = R_NilValue,
  const std::string& metric = "footrule",
  const int& leap_size = 1,
  Rcpp::Nullable<arma::cube> aug_init = R_NilValue,
  int num_obs = 0
) {

  int num_new_obs = new_rankings.n_cols;
  int n_users = rankings.n_cols;
  int n_items = rankings.n_rows;
  vec obs_freq = ones(n_users);
  vec aug_prob = ones(n_particles);

  mat rho_samples(n_items, n_particles);
  vec alpha_samples = zeros(n_particles);
  double effective_sample_size;

  cube aug_rankings;
  umat missing_indicator;
  if(type == "partial"){

    rankings.replace(datum::nan, 0);
    missing_indicator = conv_to<umat>::from(rankings);
    missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
    aug_rankings = zeros(n_items, n_users, n_particles);

    for(int i{}; i < n_particles; i++) {
      aug_rankings.slice(i) = rankings;
      initialize_missing_ranks(aug_rankings.slice(i), missing_indicator);
    }

    if(aug_init.isNotNull()) {
      aug_rankings(span::all, span(0, num_obs - 1), span::all) = Rcpp::as<cube>(aug_init);
    }
  }

  num_obs += num_new_obs;
  mat new_observed_rankings, all_observed_rankings;
  if(type == "complete"){
    new_observed_rankings = new_rankings;
    all_observed_rankings = rankings;
  }

  rho_samples = rho_init;
  alpha_samples = alpha_init;
  vec log_inc_wgt(n_particles, fill::zeros);

  if(type == "partial"){
    smc_mallows_new_users_augment_partial(
      aug_rankings, aug_prob, rho_samples, alpha_samples, num_obs, num_new_obs,
      aug_method, missing_indicator, metric);
  }


  vec norm_wgt;
  smc_mallows_new_users_reweight(
    log_inc_wgt, effective_sample_size, norm_wgt, aug_rankings, new_observed_rankings, rho_samples,
    alpha_samples, logz_estimate, cardinalities, num_obs, num_new_obs, aug_prob,
    type != "complete", metric);

  smc_mallows_new_users_resample(
    rho_samples, alpha_samples, aug_rankings, norm_wgt, num_obs,
    type != "complete");

  for (int ii = 0; ii < n_particles; ++ii) {

    for (int kk = 0; kk < mcmc_steps; ++kk) {
      if(type == "partial") {
        all_observed_rankings = aug_rankings.slice(ii);
      }
      rho_samples.col(ii) = make_new_rho(rho_samples.col(ii), all_observed_rankings,
                      alpha_samples(ii), leap_size, metric, obs_freq);


      alpha_samples(ii) = update_alpha(alpha_samples(ii), all_observed_rankings,
                    obs_freq, rho_samples.col(ii), alpha_prop_sd, metric,
                    lambda, cardinalities, logz_estimate);

    }

    if(type == "partial") {

      for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
        double log_hastings_correction = 0;

        uvec unranked_items = shuffle(find(missing_indicator.col(jj) == 1));
        if(aug_method == "uniform") {
          aug_rankings(span::all, span(jj), span(ii)) = make_new_augmentation(
            aug_rankings(span::all, span(jj), span(ii)),
            missing_indicator.col(jj),
            alpha_samples(ii),
            rho_samples.col(ii),
            metric
          );
        } else {

          Rcpp::List pprop = make_pseudo_proposal(
            unranked_items, aug_rankings(span::all, span(jj), span(ii)),
            alpha_samples(ii), rho_samples.col(ii), metric, true
          );

          Rcpp::List bprop = make_pseudo_proposal(
            unranked_items, aug_rankings(span::all, span(jj), span(ii)),
            alpha_samples(ii), rho_samples.col(ii), metric, false);
          double bprob = bprop["probability"];

          vec ar = pprop["proposal"];
          double prob = pprop["probability"];

          log_hastings_correction = -std::log(prob) + std::log(bprob);

          double ratio = -alpha_samples(ii) / n_items * (
            get_rank_distance(ar, rho_samples.col(ii), metric) -
              get_rank_distance(aug_rankings(span::all, span(jj), span(ii)),
                                rho_samples.col(ii), metric)
          ) + log_hastings_correction;


          double u = std::log(randu<double>());

          if(ratio > u){
            aug_rankings(span::all, span(jj), span(ii)) = ar;
          }
        }


      }
    }



  }

  // return the history of the particles and their values
  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = rho_samples,
    Rcpp::Named("alpha_samples") = alpha_samples,
    Rcpp::Named("augmented_rankings") = aug_rankings,
    Rcpp::Named("ESS") = effective_sample_size
  );

  return particle_history;
}
