#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "sample.h"
#include "smc.h"
#include "partitionfuns.h"
#include "misc.h"
#include "smc_mallows_new_users.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users(
  const arma::mat& rankings,
  const arma::mat& new_rankings,
  arma::mat rho_init,
  arma::vec alpha_init,
  const std::string& type,
  const int& n_particles,
  const int& mcmc_steps,
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const double alpha = 0,
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

  mat rho_samples(n_items, n_particles);
  vec alpha_samples = zeros(n_particles);
  double effective_sample_size;

  cube aug_rankings; // no. users by items by particles
  if(type == "partial" || type == "partial_alpha_fixed"){
    aug_rankings = zeros(n_users, n_items, n_particles);
    if(aug_init.isNotNull()) {
      aug_rankings(span(0, num_obs - 1), span::all, span::all) = Rcpp::as<cube>(aug_init);
    }
  }

  num_obs += num_new_obs;

  /* ====================================================== */
  /* New Information                                        */
  /* ====================================================== */
  // create two ranking dataset to use for the reweight and move stages of the
  // algorithm

  mat new_observed_rankings, all_observed_rankings;
  if(type == "complete"){
    new_observed_rankings = new_rankings;
    all_observed_rankings = rankings;
  }

  // propagate particles onto the next time step
  rho_samples = rho_init;
  alpha_samples = alpha_init;
  vec log_inc_wgt(n_particles, fill::zeros);

  /* ====================================================== */
  /* Augment partial rankings                               */
  /* ====================================================== */

  vec aug_prob = ones(n_particles);
  if(type == "partial" || type == "partial_alpha_fixed"){
    smc_mallows_new_users_augment_partial(
      aug_rankings, aug_prob, rho_samples, alpha_samples, num_obs, num_new_obs,
      rankings, aug_method, alpha, type != "partial_alpha_fixed", metric);
  }



  /* ====================================================== */
  /* Re-weight                                              */
  /* ====================================================== */

  // calculate incremental weight for each particle, based on
  // new observed rankings

  vec norm_wgt;
  smc_mallows_new_users_reweight(
    log_inc_wgt, effective_sample_size, norm_wgt, aug_rankings, new_observed_rankings, rho_samples,
    alpha, alpha_samples, logz_estimate, cardinalities, num_obs, num_new_obs, ones(n_particles),
    type != "partial_alpha_fixed", type != "complete", metric);



  /* ====================================================== */
  /* Resample                                               */
  /* ====================================================== */


  smc_mallows_new_users_resample(
    rho_samples, alpha_samples, aug_rankings, norm_wgt, num_obs,
    type != "partial_alpha_fixed", type != "complete");


  /* ====================================================== */
  /* Move step                                              */
  /* ====================================================== */
  for (int ii = 0; ii < n_particles; ++ii) {
    if(type == "complete"){
      for (int kk = 0; kk < mcmc_steps; ++kk) {
        rho_samples.col(ii) = make_new_rho(rho_samples.col(ii), all_observed_rankings,
                        alpha_samples(ii), leap_size, metric, obs_freq);

        alpha_samples(ii) = metropolis_hastings_alpha(
          alpha_samples(ii), n_items, all_observed_rankings, rho_samples.col(ii), logz_estimate,
          cardinalities, metric, alpha_prop_sd, lambda
        );
      }
    } else if(type == "partial" || type == "partial_alpha_fixed"){
        for (int kk = 0; kk < mcmc_steps; ++kk) {
          double as = (type == "partial" ? alpha_samples(ii) : alpha);
          all_observed_rankings = aug_rankings(span(0, num_obs - 1), span::all, span(ii));

          // move each particle containing sample of rho and alpha by using
          // the MCMC kernels
          rho_samples.col(ii) =
            make_new_rho(
              rho_samples.col(ii), all_observed_rankings, as, leap_size, metric, obs_freq
            );

          if(type == "partial"){
            alpha_samples(ii) = metropolis_hastings_alpha(
              as, n_items, all_observed_rankings, rho_samples.col(ii), logz_estimate,
              cardinalities, metric, alpha_prop_sd, lambda
            );
          }
          for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
            rowvec ar;
            ar = aug_rankings(span(jj), span::all, span(ii));
            vec mh_aug_result;
            mh_aug_result = metropolis_hastings_aug_ranking(
              as, rho_samples.col(ii), n_items, rankings.row(jj).t(), ar.t(),
              is_pseudo(aug_method, metric), metric
            );
            aug_rankings(span(jj), span::all, span(ii)) = mh_aug_result;
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
