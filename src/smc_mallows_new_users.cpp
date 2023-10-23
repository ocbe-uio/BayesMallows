#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "sample.h"
#include "smc.h"
#include "partitionfuns.h"
#include "misc.h"
#include "smc_mallows_new_users.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-Mallows New Users
//' @description Function to perform resample-move SMC algorithm where we
//' receive new users with complete rankings at each time step. See Chapter 4
//' of \insertCite{steinSequentialInferenceMallows2023}{BayesMallows}
//'
//' @param rankings Matrix containing the full set of observed rankings of size
//' n_assessors by n_items
//' @param type One of \code{"complete"}, \code{"partial"}, or
//' \code{"partial_alpha_fixed"}.
//' @param n_items Integer is the number of items in a ranking
//' @param n_users number of users
//' @param metric A character string specifying the distance metric to use
//' in the Bayesian Mallows Model. Available options are \code{"footrule"},
//' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//' \code{"ulam"}.
//' @param leap_size leap_size Integer specifying the step size of the
//' leap-and-shift proposal distribution
//' @param N Integer specifying the number of particles
//' @param Time Integer specifying the number of time steps in the SMC algorithm
//' @param logz_estimate Estimate of the partition function, computed with
//' \code{\link{estimate_partition_function}}.
//' @param cardinalities Cardinalities for exact evaluation of partition function,
//' returned from \code{\link{prepare_partition_function}}.
//' @param mcmc_kernel_app Integer value for the number of applications we
//' apply the MCMC move kernel
//' @param num_new_obs Integer value for the number of new observations
//' (complete rankings) for each time step
//' @param alpha_prop_sd Numeric value specifying the standard deviation of the
//'   lognormal proposal distribution used for \eqn{\alpha} in the
//'   Metropolis-Hastings algorithm. Defaults to \code{0.1}.
//' @param lambda Strictly positive numeric value specifying the rate parameter
//'   of the truncated exponential prior distribution of \eqn{\alpha}. Defaults
//'   to \code{0.1}. When \code{n_cluster > 1}, each mixture component
//'   \eqn{\alpha_{c}} has the same prior distribution.
//' @param alpha_max Maximum value of \code{alpha} in the truncated exponential
//'   prior distribution.
//' @param alpha A numeric value of the scale parameter which is known and fixed.
//' @param aug_method A character string specifying the approach for filling
//' in the missing data, options are "pseudolikelihood" or "random".
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//' @param rho_init Initial value of \code{rho}.
//'
//' @return a set of particles each containing a value of rho and alpha
//'
//' @noRd
//'
//' @example inst/examples/smc_mallows_new_users_complete_example.R
//'
//' @family modeling
//'
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users_cpp(
  const arma::mat& rankings,
  const std::string& type,
  const int& n_items,
  const int& n_users,
  const int& N,
  int Time,
  const int& mcmc_kernel_app,
  const int& num_new_obs,
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const double alpha_max = 1e6,
  const double alpha = 0,
  const std::string& aug_method = "random",
  const Rcpp::Nullable<arma::vec>& logz_estimate = R_NilValue,
  const Rcpp::Nullable<arma::vec>& cardinalities = R_NilValue,
  const bool verbose = false,
  const std::string& metric = "footrule",
  const int& leap_size = 1,
  Rcpp::Nullable<arma::mat> rho_init = R_NilValue,
  Rcpp::Nullable<arma::vec> alpha_init = R_NilValue
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */

  /* generate rho samples using uniform prior ------------- */
  cube rho_samples(N, n_items, Time + 1);
  rho_samples.slice(0) = initialize_rho(n_items, N, rho_init).t();

  /* generate alpha samples using exponential prior ------- */
  mat alpha_samples;
  if(type != "partial_alpha_fixed") {
    alpha_samples = zeros(N, Time + 1);
    alpha_samples.col(0) = initialize_alpha(N, alpha_init);
  }

  /* generate vector to store ESS */
  rowvec ESS_vec(Time);

  // this is to store the augmentations of the observed rankings for each particle
  cube aug_rankings; // no. users by items by particles
  if(type == "partial" || type == "partial_alpha_fixed"){
    aug_rankings = zeros(n_users, n_items, N);
  }

  /* ====================================================== */
  /* New user situation                                     */
  /* ====================================================== */
  int num_obs = 0;

  for (int tt{}; tt < Time; ++tt) {
    if (verbose) REprintf("observe %i out of %i \n", tt + 1, Time);

    // keep tally of how many ranking observations we have so far
    num_obs += num_new_obs;

    /* ====================================================== */
    /* New Information                                        */
    /* ====================================================== */
    // create two ranking dataset to use for the reweight and move stages of the
    // algorithm
    int row_start;
    mat new_observed_rankings, all_observed_rankings;
    if(type == "complete"){
      row_start = num_obs - num_new_obs;
      new_observed_rankings = rankings.submat(row_start, 0, num_obs - 1, rankings.n_cols - 1);
      all_observed_rankings = rankings.submat(0, 0, num_obs - 1, rankings.n_cols - 1);
    }

    // propagate particles onto the next time step
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    if(type != "partial_alpha_fixed") alpha_samples.col(tt + 1) = alpha_samples.col(tt);
    vec log_inc_wgt(N, fill::zeros);

    /* ====================================================== */
    /* Augment partial rankings                               */
    /* ====================================================== */

    vec aug_prob = ones(N);
    if(type == "partial" || type == "partial_alpha_fixed"){
      smc_mallows_new_users_augment_partial(
        aug_rankings, aug_prob, rho_samples, alpha_samples, num_obs, num_new_obs,
        rankings, aug_method, tt, alpha, type != "partial_alpha_fixed", metric);
    }

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    // calculate incremental weight for each particle, based on
    // new observed rankings

    vec norm_wgt;
    smc_mallows_new_users_reweight(
      log_inc_wgt, ESS_vec, norm_wgt, aug_rankings, new_observed_rankings, rho_samples,
      alpha, alpha_samples, tt, logz_estimate, cardinalities, num_obs, num_new_obs, ones(N),
      type != "partial_alpha_fixed", type != "complete", metric);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */

    smc_mallows_new_users_resample(
      rho_samples, alpha_samples, aug_rankings, norm_wgt, tt, num_obs,
      type != "partial_alpha_fixed", type != "complete");

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (int ii = 0; ii < N; ++ii) {
      if(type == "complete"){
        for (int kk = 0; kk < mcmc_kernel_app; ++kk) {
          // move each particle containing sample of rho and alpha by using
          // the MCMC kernels
          const double& as = alpha_samples(ii, tt + 1);
          const rowvec& rs = \
            rho_samples(span(ii), span::all, span(tt + 1));
          rho_samples(span(ii), span::all, span(tt + 1)) =
            metropolis_hastings_rho(
              as, n_items, all_observed_rankings, rs.t(), metric, leap_size
            );
          alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(
            as, n_items, all_observed_rankings, rs.t(), logz_estimate,
            cardinalities, metric, alpha_prop_sd, alpha_max, lambda
          );
        }
      } else if(type == "partial" || type == "partial_alpha_fixed"){
          for (int kk = 0; kk < mcmc_kernel_app; ++kk) {
            double as = (type == "partial" ? alpha_samples(ii, tt + 1) : alpha);
            mat all_observed_rankings;
            all_observed_rankings = aug_rankings(span(0, num_obs - 1), span::all, span(ii));
            mat rs_slice = rho_samples.slice(tt + 1);
            rowvec rs = rs_slice.row(ii);
            // move each particle containing sample of rho and alpha by using
            // the MCMC kernels
            rho_samples(span(ii), span::all, span(tt + 1)) =
              metropolis_hastings_rho(
                as, n_items, all_observed_rankings, rs.t(), metric, leap_size
              );
            if(type == "partial"){
              alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(
                as, n_items, all_observed_rankings, rs.t(), logz_estimate,
                cardinalities, metric, alpha_prop_sd, alpha_max, lambda
              );
            }
            for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
              rowvec ar;
              ar = aug_rankings(span(jj), span::all, span(ii));
              vec mh_aug_result;
              mh_aug_result = metropolis_hastings_aug_ranking(
                as, rs.t(), n_items, rankings.row(jj).t(), ar.t(),
                is_pseudo(aug_method, metric), metric
              );
              aug_rankings(span(jj), span::all, span(ii)) = mh_aug_result;
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
    Rcpp::Named("ESS") = ESS_vec
  );
  particle_history.attr("class") = "SMCMallows";
  return particle_history;
}
