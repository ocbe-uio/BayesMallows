#include "RcppArmadillo.h"
#include "misc.h"
#include "smc.h"
#include "partitionfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-mallows new users partial (alpha fixed)
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
//' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
//' @param alpha A numeric value of the scale parameter which is known and fixed
//' @return a set of particles each containing a value of rho and alpha
//' @export
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_users_partial_alpha_fixed(
  const arma::mat& R_obs,
  const unsigned int& n_items,
  const std::string metric,
  const int& leap_size,
  const unsigned int& N,
  const unsigned int Time,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const int& mcmc_kernel_app,
  const unsigned int& num_new_obs,
  const std::string& aug_method,
  const double alpha
) {

  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */
  const int& n_users = R_obs.n_rows; // this is total- number of users

  // generate rho samples using uniform prior
  arma::cube rho_samples(N, n_items, Time + 1, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) {
    const arma::uvec items_sample = arma::randperm(n_items) + 1;
    for (arma::uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  // this is to store the augmentations of the observed rankings for each particle
  arma::cube aug_rankings(n_users, n_items, N, arma::fill::zeros); // no. users by items by particles

  /* ====================================================== */
  /* New user situation                                     */
  /* ====================================================== */
  unsigned int num_obs = 0;

  for (arma::uword tt = 0; tt < Time; ++tt) {

    /* ====================================================== */
    /* New Information                                        */
    /* ====================================================== */
    // keep tally of how many ranking observations we have so far
    num_obs += num_new_obs;

    // create two ranking dataset to use for the reweight and move stages of the algorithm
    // Note:
    // new_observed_rankings = R_obs[((num_obs-num_new_obs+1):num_obs),]
    // all_observed_rankings = R_obs[(1:num_obs),]

    // propagate particles onto the next time step
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);

    // calculate incremental weight and augmentation prob for each particle,
    // based on new observed rankings
    arma::vec log_inc_wgt(N, arma::fill::zeros);

    /* ====================================================== */
    /* Augment partial rankings                               */
    /* ====================================================== */

    const arma::ivec ranks = Rcpp::seq(1, n_items);
    arma::vec aug_prob = Rcpp::rep(1.0, N);

    for (arma::uword ii = 0; ii < N; ++ii) {
      for (arma::uword jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
        arma::vec partial_ranking = R_obs.row(jj).t();

        // find items missing from original observed ranking
        const arma::uvec& unranked_items = find_nonfinite(partial_ranking);

        // find ranks missing from ranking
        const Rcpp::NumericVector& missing_ranks = Rcpp_setdiff_arma(ranks, partial_ranking);

        // fill in missing ranks based on choice of augmentation method
        if (aug_method == "random") {
        // create new agumented ranking by sampling remaining ranks from set uniformly
          if (missing_ranks.length() == 1) {
            partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = Rcpp::as<arma::vec>(missing_ranks);
          } else {
            partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = Rcpp::as<arma::vec>(Rcpp::sample(missing_ranks, missing_ranks.length()));
          }

          aug_rankings(arma::span(jj), arma::span::all, arma::span(ii)) = partial_ranking;
          aug_prob(ii) = divide_by_fact(aug_prob(ii), missing_ranks.length());
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          // randomly permute the unranked items to give the order in which they will be allocated
          arma::uvec item_ordering;
          item_ordering = arma::conv_to<arma::uvec>::from(arma::shuffle(unranked_items));
          const arma::rowvec& rho_s = rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1));
          const Rcpp::List& proposal = calculate_forward_probability(\
            item_ordering, partial_ranking, missing_ranks, rho_s.t(),\
            alpha, n_items, metric\
          );
          const arma::vec& a_rank = proposal["aug_ranking"];
          const double& f_prob = proposal["forward_prob"];
          aug_rankings(arma::span(jj), arma::span::all, arma::span(ii)) = a_rank;
          aug_prob(ii) = aug_prob(ii) * f_prob;
        } else {
          Rcpp::stop("Combined choice of metric and aug_method is incompatible");
        }
      }
    }

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    for (arma::uword ii = 0; ii < N; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Initializing variables ------------------------------- */
      const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue;
      const arma::rowvec rho_samples_ii = \
        rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1));

      /* Calculating log_z_alpha and log_likelihood ----------- */
      const double log_z_alpha = get_partition_function(\
        n_items, alpha, cardinalities, logz_estimate, metric\
      );

      arma::mat new_observed_rankings;
      new_observed_rankings = aug_rankings(arma::span(num_obs - num_new_obs, num_obs - 1), arma::span::all, arma::span(ii));
      double log_likelihood = get_mallows_loglik(\
        alpha, rho_samples_ii.t(), n_items, new_observed_rankings, metric\
      );
      log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
    }

    /* normalise weights ------------------------------------ */
    const double maxw = arma::max(log_inc_wgt);
    const arma::vec w = arma::exp(log_inc_wgt - maxw);
    const arma::vec norm_wgt = w / arma::sum(w);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */

    /* Resample particles using multinomial resampling ------ */
    arma::uvec index = permute_with_weights(norm_wgt, N);
    arma::uvec tt_vec;
    tt_vec = tt;

    /* Replacing tt + 1 slice on rho_samples ---------------- */
    arma::mat rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
    rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
    rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

    /* Replacing tt + 1 column on alpha_samples ------------- */
    const arma::cube& aug_rankings_index = aug_rankings.slices(index);
    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(arma::span(0, num_obs - 1), arma::span::all, arma::span::all);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (arma::uword ii = 0; ii < N; ++ii) {
      arma::mat all_observed_rankings;
      all_observed_rankings = aug_rankings(arma::span(0, num_obs - 1), arma::span::all, arma::span(ii));
      const arma::mat& rs_slice = rho_samples.slice(tt + 1);
      const arma::rowvec& rs = rs_slice.row(ii);
      // move each particle containing sample of rho and alpha by using
      // the MCMC kernels
      rho_samples(arma::span(ii), arma::span::all, arma::span(tt + 1)) =\
        metropolis_hastings_rho(\
          alpha, n_items, all_observed_rankings, metric, rs.t(), leap_size\
        );
      for (arma::uword jj = 0; jj < num_obs; ++jj) {
        arma::rowvec ar;
        ar = aug_rankings(arma::span(jj), arma::span::all, arma::span(ii));
        arma::vec mh_aug_result;
        if (aug_method == "random") {
          mh_aug_result = metropolis_hastings_aug_ranking(\
          alpha, rs.t(), n_items, R_obs.row(jj).t(), ar.t(), metric\
        );
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          mh_aug_result = metropolis_hastings_aug_ranking_pseudo(
            alpha, rs.t(), n_items, R_obs.row(jj).t(), ar.t(), metric\
          );
        }
        aug_rankings(arma::span(jj), arma::span::all, arma::span(ii)) = mh_aug_result;
      }
    }
  }
  // return the history of the particles and their values
  return Rcpp::List::create((Rcpp::Named("rho_samples") = rho_samples));
}
