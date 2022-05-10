#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"
#include "smc_mallows_new_users.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-Mallows new users partial
//' @description Function to perform resample-move SMC algorithm where we receive new users with complete rankings
//' at each time step
//' @param R_obs Matrix containing the full set of observed rankings of size n_assessors by n_items
//' @param n_items Integer is the number of items in a ranking
//' @param metric A character string specifying the distance metric to use in the
//' Bayesian Mallows Model. Available options are \code{"footrule"},
//' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//' \code{"ulam"}.
//' @param leap_size leap_size Integer specifying the step size of the leap-and-shift
//' proposal distribution
//' @param N Integer specifying the number of particles
//' @param Time Integer specifying the number of time steps in the SMC algorithm
//' @param logz_estimate Estimate of the partition function, computed with
//' \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
//' @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel
//' @param num_new_obs Integer value for the number of new observations (complete rankings) for each time step
//' @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha
//' @param lambda Strictly positive numeric value specifying the rate parameter
//' of the truncated exponential prior distribution of alpha.
//' @param alpha_max  Maximum value of alpha in the truncated exponential
//' prior distribution.
//' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//' @return a set of particles each containing the values of rho and alpha and the effective sample size (ESS) at each iteration of the SMC
//' algorithm as well as the set of augmented rankings at the final iteration.
//' @export
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users_partial(
  const arma::mat& R_obs,
  const int& n_items,
  const std::string& metric,
  const int& leap_size,
  const int& N,
  const int Time,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const int& mcmc_kernel_app,
  const int& num_new_obs,
  const double alpha_prop_sd,
  const double lambda,
  const double alpha_max,
  const std::string& aug_method,
  const bool verbose = false
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */
  int n_users = R_obs.n_rows; // this is total- number of users

  /* generate rho samples using uniform prior ------------- */
  cube rho_samples = initialize_rho(N, n_items, Time + 1);

  /* generate alpha samples using exponential prior ------- */
  mat alpha_samples(N, Time + 1);
  const vec& alpha_samples_0 = Rcpp::rexp(N, 1);
  alpha_samples.col(0) = alpha_samples_0;

  /* generate vector to store ESS */
  rowvec ESS_vec(Time);

  // this is to store the augmentations of the observed rankings for each particle
  cube aug_rankings(n_users, n_items, N, fill::zeros); // no. users by items by particles

  /* ====================================================== */
  /* New user situation                                     */
  /* ====================================================== */
  unsigned int num_obs = 0;

  for (uword tt = 0; tt < Time; ++tt) {
    if (verbose) REprintf("observe %i out of %i \n", tt + 1, Time);
    /* ====================================================== */
    /* New Information                                        */
    /* ====================================================== */
    // keep tally of how many ranking observations we have so far
    num_obs += num_new_obs;

    // propagate particles onto the next time step
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    alpha_samples.col(tt + 1) = alpha_samples.col(tt);

    // calculate incremental weight and augmentation prob for each particle,
    // based on new observed rankings
    vec log_inc_wgt(N, fill::zeros);

    /* ====================================================== */
    /* Augment partial rankings                               */
    /* ====================================================== */

    vec aug_prob = ones(N);
    smc_mallows_new_users_augment_partial(aug_rankings, aug_prob,
      rho_samples, alpha_samples, num_obs, num_new_obs, R_obs,
      aug_method, metric, tt, 0, true);

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    vec norm_wgt;
    smc_mallows_new_users_reweight(
      log_inc_wgt, ESS_vec, norm_wgt, aug_rankings, mat(), rho_samples,
      0, alpha_samples, tt, logz_estimate, metric, num_obs,
      num_new_obs, aug_prob, true, true);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */

    smc_mallows_new_users_resample(
      rho_samples, alpha_samples, aug_rankings, norm_wgt, tt, num_obs, true, true);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (uword ii = 0; ii < N; ++ii) {
      double as = alpha_samples(ii, tt + 1);
      mat all_observed_rankings;
      all_observed_rankings = aug_rankings(span(0, num_obs - 1), span::all, span(ii));
      mat rs_slice = rho_samples.slice(tt + 1);
      rowvec rs = rs_slice.row(ii);
      // move each particle containing sample of rho and alpha by using
      // the MCMC kernels
      rho_samples(span(ii), span::all, span(tt + 1)) =\
        metropolis_hastings_rho(\
          as, n_items, all_observed_rankings, metric, rs.t(), leap_size\
        );
      alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(\
        as, n_items, all_observed_rankings, metric, rs.t(), logz_estimate,\
        alpha_prop_sd, lambda, alpha_max\
      );
      for (uword jj = 0; jj < num_obs; ++jj) {
        rowvec ar;
        ar = aug_rankings(span(jj), span::all, span(ii));
        vec mh_aug_result;
        if (aug_method == "random") {
          mh_aug_result = metropolis_hastings_aug_ranking(as, rs.t(), n_items, R_obs.row(jj).t(), ar.t(), metric);
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          mh_aug_result = metropolis_hastings_aug_ranking_pseudo(as, rs.t(), n_items, R_obs.row(jj).t(), ar.t(), metric);
        }
        aug_rankings(span(jj), span::all, span(ii)) = mh_aug_result;
      }
    }
  }
  // return the history of the particles and their values
  return Rcpp::List::create(
    Rcpp::Named("rho_samples") = rho_samples,
    Rcpp::Named("alpha_samples") = alpha_samples,
    Rcpp::Named("augmented_rankings") = aug_rankings,
    Rcpp::Named("ESS") = ESS_vec
  );
}
