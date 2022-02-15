#include "RcppArmadillo.h"
#include "smc.h"
#include "partitionfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-Mallows New Users Complete
//' @description Function to perform resample-move SMC algorithm where we
//' receive new users with complete rankings at each time step
//'
//' @param R_obs Matrix containing the full set of observed rankings of size
//' n_assessors by n_items
//' @param n_items Integer is the number of items in a ranking
//' @param metric A character string specifying the distance metric to use
//' in the Bayesian Mallows Model. Available options are \code{"footrule"},
//' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//' \code{"ulam"}.
//' @param leap_size leap_size Integer specifying the step size of the
//' leap-and-shift proposal distribution
//' @param N Integer specifying the number of particles
//' @param Time Integer specifying the number of time steps in the SMC algorithm
//' @param logz_estimate Estimate of the partition function, computed with
//' \code{\link{estimate_partition_function}} in the BayesMallow R package
//' {estimate_partition_function}.
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
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//'
//' @return a set of particles each containing a value of rho and alpha
//'
//' @importFrom stats rexp
//' @export
//'
//' @example inst/examples/smc_mallows_new_users_complete_example.R
//'
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users_complete(
  const arma::mat& R_obs,
  const int& n_items,
  const std::string& metric,
  const int& leap_size,
  const int& N,
  int Time,
  const int& mcmc_kernel_app,
  const int& num_new_obs,
  const double alpha_prop_sd,
  const double lambda,
  const double alpha_max,
  const Rcpp::Nullable<arma::vec>& logz_estimate = R_NilValue,
  const bool verbose = false
) {

  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */
  const int n_users = R_obs.n_rows; // total number of users
  if (Time > n_users / num_new_obs) {
    Rcpp::warning(\
      "Time should not exceed n_users / num_new_obs. Recalculating."\
    );
    Time = n_users / num_new_obs;
  }

  /* generate rho samples using uniform prior ------------- */
  arma::cube rho_samples(N, n_items, (n_users + Time + 1), arma::fill::zeros);
  for (int i = 0; i < N; ++i) {
    const arma::uvec items_sample = arma::randperm(n_items) + 1;
    for (int j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  /* generate alpha samples using exponential prior ------- */
  arma::mat alpha_samples(N, (n_users + Time + 1));
  const arma::vec alpha_samples_0 = Rcpp::rexp(N, 1);
  alpha_samples.col(0) = alpha_samples_0;

  /* ====================================================== */
  /* New user situation                                     */
  /* ====================================================== */
  int num_obs = 0;

  for (arma::uword tt = 0; tt < Time; ++tt) {
    if (verbose) REprintf("observe %i out of %i \n", tt + 1, Time);

    // keep tally of how many ranking observations we have so far
    num_obs += num_new_obs;

    /* ====================================================== */
    /* New Information                                        */
    /* ====================================================== */
    // create two ranking dataset to use for the reweight and move stages of the
    // algorithm
    const int row_start = num_obs - num_new_obs;
    arma::mat new_observed_rankings(num_obs, R_obs.n_cols);
    arma::mat all_observed_rankings;
    new_observed_rankings = R_obs.submat(row_start, 0, num_obs - 1, R_obs.n_cols - 1);
    all_observed_rankings = R_obs.submat(0, 0, num_obs - 1, R_obs.n_cols - 1);

    // propagate particles onto the next time step
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    alpha_samples.col(tt + 1) = alpha_samples.col(tt);

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    // calculate incremental weight for each particle, based on
    // new observed rankings
    arma::vec log_inc_wgt(N, arma::fill::zeros);

    for (int ii = 0; ii < N; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Initializing variables ------------------------------- */
      const Rcpp::Nullable<arma::vec>& cardinalities = R_NilValue;
      const double& alpha_samples_ii = alpha_samples(ii, tt + 1);
      const arma::rowvec rho_samples_ii = \
        rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1));

      /* Calculating log_z_alpha and log_likelihood ----------- */
      double log_z_alpha, log_likelihood;
      log_z_alpha = get_partition_function(\
        n_items, alpha_samples_ii, cardinalities, logz_estimate, metric\
      );
      log_likelihood = get_mallows_loglik(\
        alpha_samples_ii, rho_samples_ii.t(), n_items, new_observed_rankings,\
        metric\
      );
      log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha;
    }

    /* normalise weights ------------------------------------ */
    const double& maxw = arma::max(log_inc_wgt);
    const arma::vec& w = arma::exp(log_inc_wgt - maxw);
    const arma::vec norm_wgt = w / arma::sum(w);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */

    /* Resample particles using multinomial resampling ------ */
    arma::uvec index = permutate_with_weights(norm_wgt, N);
    arma::uvec tt_vec;
    tt_vec = tt;

    /* Replacing tt + 1 slice on rho_samples ---------------- */
    arma::mat& rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
    rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
    rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

    /* Replacing tt + 1 column on alpha_samples ------------- */
    alpha_samples.col(tt + 1) =  alpha_samples.submat(index, tt_vec + 1);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (int ii = 0; ii < N; ++ii) {
      for (int kk = 0; kk < mcmc_kernel_app; ++kk) {
        // move each particle containing sample of rho and alpha by using
        // the MCMC kernels
        const double& as = alpha_samples(ii, tt + 1);
        const arma::rowvec& rs = \
          rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1));
        rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1)) =\
          metropolis_hastings_rho(\
            as, n_items, all_observed_rankings, metric, rs.t(), leap_size\
          );
        alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(\
          as, n_items, all_observed_rankings, metric, rs.t(), logz_estimate,\
          alpha_prop_sd, lambda, alpha_max\
        );
      }
    }
  }
  // return the history of the particles and their values
  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = rho_samples,
    Rcpp::Named("alpha_samples") = alpha_samples
  );
  particle_history.attr("class") = "SMCMallows"; // TODO: add List
  return particle_history;
}
