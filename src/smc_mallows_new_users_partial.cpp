#include "RcppArmadillo.h"
#include "misc.h"
#include "smc.h"
#include "partitionfuns.h"

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
//' @return a set of particles each containing a value of rho and alpha
//' @export
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users_partial(
  const arma::mat& R_obs,
  const unsigned int& n_items,
  const std::string& metric,
  const int& leap_size,
  const unsigned int& N,
  const unsigned int Time,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const int& mcmc_kernel_app,
  const unsigned int& num_new_obs,
  const double alpha_prop_sd,
  const double lambda,
  const double alpha_max,
  const std::string& aug_method,
  const bool verbose = false
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */
  const int& n_users = R_obs.n_rows; // this is total- number of users

  /* generate rho samples using uniform prior ------------- */

  cube rho_samples(N, n_items, Time + 1, fill::zeros);
  for (uword i = 0; i < N; ++i) {
    uvec items_sample = randperm(n_items) + 1;
    for (uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  /* generate alpha samples using exponential prior ------- */
  mat alpha_samples(N, Time + 1);
  const vec& alpha_samples_0 = Rcpp::rexp(N, 1);
  alpha_samples.col(0) = alpha_samples_0;

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

    const ivec ranks = Rcpp::seq(1, n_items);
    vec aug_prob = Rcpp::rep(1.0, N);

    for (uword ii = 0; ii < N; ++ii) {
      for (uword jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
        vec partial_ranking = R_obs.row(jj).t();

        // find items missing from original observed ranking
        uvec unranked_items = find_nonfinite(partial_ranking);

        // find ranks missing from ranking
        const Rcpp::NumericVector missing_ranks = Rcpp::sort_unique(Rcpp_setdiff_arma(ranks, partial_ranking));

        // fill in missing ranks based on choice of augmentation method
        if (aug_method == "random") {

          // create new agumented ranking by sampling remaining ranks from set uniformly
          if (missing_ranks.length() == 1) {
            partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = Rcpp::as<vec>(missing_ranks);
          } else {
            partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = Rcpp::as<vec>(Rcpp::sample(missing_ranks, missing_ranks.length()));
          }

          aug_rankings(span(jj), span::all, span(ii)) = partial_ranking;
          aug_prob(ii) = divide_by_fact(aug_prob(ii), missing_ranks.length());

        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {

          // randomly permute the unranked items to give the order in which they will be allocated
          uvec item_ordering;
          item_ordering = arma::conv_to<uvec>::from(arma::shuffle(unranked_items));

          const rowvec rho_s = rho_samples(span(ii), span::all, span(tt + 1));
          const Rcpp::List proposal = calculate_forward_probability(\
            item_ordering, partial_ranking, missing_ranks, rho_s.t(),\
            alpha_samples(ii, tt + 1), n_items, metric\
          );
          const vec& a_rank = proposal["aug_ranking"];
          const double& f_prob = proposal["forward_prob"];
          aug_rankings(span(jj), span::all, span(ii)) = a_rank;
          aug_prob(ii) = aug_prob(ii) * f_prob;
        } else {
          Rcpp::stop("Combined choice of metric and aug_method is incompatible");
        }
      }
    }

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    for (uword ii = 0; ii < N; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Initializing variables ------------------------------- */
      const Rcpp::Nullable<vec>& cardinalities = R_NilValue;
      double alpha_samples_ii = alpha_samples(ii, tt + 1);
      rowvec rho_samples_ii = \

        rho_samples(span(ii), span::all, span(tt + 1));

      /* Calculating log_z_alpha and log_likelihood ----------- */
      double log_z_alpha, log_likelihood;
      log_z_alpha = get_partition_function(\
        n_items, alpha_samples_ii, cardinalities, logz_estimate, metric\
      );

      mat new_observed_rankings;
      new_observed_rankings = aug_rankings(span(num_obs - num_new_obs, num_obs - 1), span::all, span(ii));
      log_likelihood = get_mallows_loglik(\
        alpha_samples_ii, rho_samples_ii.t(), n_items, new_observed_rankings,\
        metric\
      );
      log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha;
    }

    /* normalise weights ------------------------------------ */
    double maxw = arma::max(log_inc_wgt);
    vec w = arma::exp(log_inc_wgt - maxw);
    vec norm_wgt = w / arma::sum(w);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */

    /* Resample particles using multinomial resampling ------ */
    uvec index = permute_with_weights(norm_wgt, N);
    uvec tt_vec;
    tt_vec = tt;

    /* Replacing tt + 1 slice on rho_samples ---------------- */
    mat rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
    rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
    rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

    /* Replacing tt + 1 column on alpha_samples ------------- */
    alpha_samples.col(tt + 1) =  alpha_samples.submat(index, tt_vec + 1);
    cube aug_rankings_index = aug_rankings.slices(index);

    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(span(0, num_obs - 1), span::all, span::all);

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
    Rcpp::Named("alpha_samples") = alpha_samples
  );
}
