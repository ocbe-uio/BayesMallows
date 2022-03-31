#include "RcppArmadillo.h"
#include "misc.h"
#include "smc.h"
#include "partitionfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-Mallows new users rank
//' @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user
//' at each time step. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.
//' @param n_items Integer is the number of items in a ranking
//' @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time time steps
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
//' @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha
//' @param lambda Strictly positive numeric value specifying the rate parameter
//' of the truncated exponential prior distribution of alpha.
//' @param alpha_max  Maximum value of alpha in the truncated exponential
//' prior distribution.
//' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//' @return a 3d matrix containing the samples of rho and alpha from the SMC algorithm
//' @export
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_item_rank(
  const unsigned int& n_items,
  arma::cube& R_obs,
  const std::string& metric,
  const int& leap_size,
  const unsigned int& N,
  const unsigned int Time,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const int& mcmc_kernel_app,
  const double alpha_prop_sd,
  const double lambda,
  const double alpha_max,
  const std::string& aug_method,
  const bool verbose = false
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */

  // Generate N initial samples of rho using the uniform prior
  arma::cube rho_samples(N, n_items, Time, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) {
    const arma::uvec items_sample = arma::randperm(n_items) + 1;
    for (arma::uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  /* generate alpha samples using exponential prior ------- */
  arma::mat alpha_samples(N, Time);
  const arma::vec alpha_samples_0 = Rcpp::rexp(N, 1);
  alpha_samples.col(0) = alpha_samples_0;

  /* ====================================================== */
  /* Augment Rankings                                       */
  /* ====================================================== */
  const unsigned int& num_ranks = R_obs.n_rows;

  // each particle has its own set of augmented rankings
  arma::cube aug_rankings(num_ranks, n_items, N, arma::fill::zeros);
  arma::cube prev_aug_rankings(num_ranks, n_items, N, arma::fill::zeros);

  // augment incomplete ranks to initialise
  const arma::ivec ranks = Rcpp::seq(1, n_items);

  // total correction prob
  arma::vec total_correction_prob = Rcpp::rep(1.0, N);

  // iterate through each observed ranking and create new "corrected" augmented rankings
  for (arma::uword ii = 0; ii < N; ++ii) {
    // set t-1 generation to old as we sample for t new
    prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

    // make the correction
    for (arma::uword jj = 0; jj < num_ranks; ++jj) {
      // fill in missing ranks based on choice of augmentation method
      arma::vec R_obs_slice_0_row_jj = R_obs.slice(0).row(jj).t();
      if (aug_method == "random") {
        // find elements missing from original observed ranking
        const Rcpp::NumericVector remaining_set = Rcpp_setdiff_arma(ranks, R_obs_slice_0_row_jj);

        // create new agumented ranking by sampling remaining ranks from set uniformly
        arma::vec rset;
        const int& remaining_set_length = remaining_set.length();
        if (remaining_set_length == 1) {
          rset = Rcpp::as<arma::vec>(remaining_set);
        } else {
          rset = Rcpp::as<arma::vec>(Rcpp::sample(remaining_set, remaining_set.length()));
        }
        arma::vec partial_ranking = R_obs_slice_0_row_jj;
        partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = rset;

        aug_rankings.slice(ii).row(jj) = partial_ranking.t();
        total_correction_prob(ii) = divide_by_fact(total_correction_prob(ii), remaining_set_length);
      } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
        // find items missing from original observed ranking
        const arma::uvec& unranked_items = arma::find_nonfinite(R_obs_slice_0_row_jj);

        // find unallocated ranks from original observed ranking
        const Rcpp::NumericVector& remaining_set = Rcpp_setdiff_arma(ranks, R_obs_slice_0_row_jj);

        // randomly permute the unranked items to give the order in which they will be allocated
        arma::uvec item_ordering;
        item_ordering = arma::conv_to<arma::uvec>::from(arma::shuffle(unranked_items));
        const Rcpp::List proposal = calculate_forward_probability(\
          item_ordering, R_obs_slice_0_row_jj, remaining_set, rho_samples.slice(0).row(ii).t(),\
          alpha_samples(ii, 0), n_items, metric\
        );
        const arma::vec& a_rank = proposal["aug_ranking"];
        const double& f_prob = proposal["forward_prob"];
        aug_rankings(arma::span(jj), arma::span::all, arma::span(ii)) = a_rank;
        total_correction_prob(ii) *= f_prob;
      } else {
        Rcpp::stop(\
          "Combined choice of metric and aug_method is incompatible. ",
          "The value is TRUE, so the script must end here"\
        );
      }
    }
  }

  /* ====================================================== */
  /* Re-weight                                              */
  /* ====================================================== */

  // incremental weight for each particle, based on new observed rankings
  arma::vec log_inc_wgt(N, arma::fill::zeros);

  for (arma::uword ii = 0; ii < N; ++ii) {
    // evaluate the log estimate of the partition function for a particular
    // value of alpha

    /* Initializing variables ------------------------------- */
    const Rcpp::Nullable<arma::vec>& cardinalities = R_NilValue;

    /* Calculating log_z_alpha and log_likelihood ----------- */
    const double& log_z_alpha = get_partition_function(\
      n_items, alpha_samples(ii, 0), cardinalities, logz_estimate, metric\
    );
    const double& log_likelihood = get_mallows_loglik(\
      alpha_samples(ii, 0), rho_samples.slice(0).row(ii).t(), n_items,\
      aug_rankings.slice(ii), metric\
    );
    const double& log_tcp = std::log(total_correction_prob(ii));
    log_inc_wgt(ii) = log_likelihood - num_ranks * log_z_alpha - log_tcp;
  }

    /* normalise weights ------------------------------------ */
    const double& maxw = arma::max(log_inc_wgt);
    const arma::vec& w = arma::exp(log_inc_wgt - maxw);
    const arma::vec& norm_wgt = w / arma::sum(w);

  /* ====================================================== */
  /* Resample                                               */
  /* ====================================================== */
  /* Resample particles using multinomial resampling ------ */
  arma::uvec index = permutate_with_weights(norm_wgt, N);
  rho_samples.slice(0) = rho_samples.slice(0).rows(index);
  const arma::vec& asc = alpha_samples.col(0);
  alpha_samples.col(0) = asc.elem(index);
  aug_rankings = aug_rankings.slices(index);

  /* ====================================================== */
  /* Move step                                              */
  /* ====================================================== */
  for (arma::uword ii = 0; ii < N; ++ii) {
    rho_samples.slice(0).row(ii) = metropolis_hastings_rho(\
      alpha_samples(ii, 0), n_items, aug_rankings.slice(ii), metric,\
      rho_samples.slice(0).row(ii).t(), leap_size\
    ).t();
    alpha_samples(ii, 0) = metropolis_hastings_alpha(\
      alpha_samples(ii, 0), n_items, aug_rankings.slice(ii), metric,\
        rho_samples.slice(0).row(ii).t(), logz_estimate,\
        alpha_prop_sd, lambda, alpha_max\
    );
    for (arma::uword jj = 0; jj < num_ranks; ++jj) {
      arma::vec mh_aug_result;
      if (aug_method == "random") {
        mh_aug_result = metropolis_hastings_aug_ranking(\
          alpha_samples(ii, 0), rho_samples.slice(0).row(ii).t(), n_items,\
          R_obs.slice(0).row(jj).t(), aug_rankings.slice(ii).row(jj).t(), metric\
        );
      } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
        mh_aug_result = metropolis_hastings_aug_ranking_pseudo(\
          alpha_samples(ii, 0), rho_samples.slice(0).row(ii).t(), n_items,\
          R_obs.slice(0).row(jj).t(), aug_rankings.slice(ii).row(jj).t(), metric\
        );
      }
      aug_rankings.slice(ii).row(jj) = mh_aug_result.t();
    }
  }

  /* ====================================================== */
  /* Loop for t=1,...,Time                                  */
  /* ====================================================== */

  for (arma::uword tt = 0; tt < Time - 1; ++tt) {
    if (verbose) REprintf("iteration %i out of %i \n", tt + 1, Time - 1);

    /* New Information -------------------------------------- */
    // new observed item ranks from each user, need to update augmented rankings
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    alpha_samples.col(tt + 1) = alpha_samples.col(tt);

    // total correction prob
    arma::vec particle_correction_prob = Rcpp::rep(1.0, N);

    // iterate through each observed ranking and create new "corrected"
    // augmented rankings

    for (arma::uword ii = 0; ii < N; ++ii) {
      // set t-1 generation to old as we sample for t new
      prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

      // make the correction
      for (arma::uword jj = 0; jj < num_ranks; ++jj) {
        if (aug_method == "random") {
          const Rcpp::List check_correction = correction_kernel(\
            R_obs.slice(tt + 1).row(jj).t(), aug_rankings.slice(ii).row(jj).t(),\
            n_items
          );
          const arma::vec c_rank = check_correction["ranking"];
          const double c_prob = check_correction["correction_prob"];
          aug_rankings.slice(ii).row(jj) = c_rank.t();
          particle_correction_prob(ii) *= c_prob;
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          const Rcpp::List check_correction = correction_kernel_pseudo(\
            aug_rankings.slice(ii).row(jj).t(), R_obs.slice(tt + 1).row(jj).t(),\
            rho_samples.slice(tt + 1).row(ii).t(), alpha_samples(ii, tt + 1),\
            n_items, metric
          );
          const arma::vec c_rank = check_correction["ranking"];
          const double c_prob = check_correction["correction_prob"];
          aug_rankings.slice(ii).row(jj) = c_rank.t();
          // # these probs are in real scale
          particle_correction_prob(ii) *= c_prob;
        } else {
          Rcpp::stop("Combined choice of metric and aug_method is incompatible");
        }
      }
    }

    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    // incremental weight for each particle, based on new observed rankings
    arma::vec log_inc_wgt(N, arma::fill::zeros);
    for (arma::uword ii = 0; ii < N; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Calculating log_z_alpha and log_likelihood ----------- */
      double loglik_1 = get_mallows_loglik(\
        alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        aug_rankings.slice(ii), metric\
      );
      double loglik_2 = get_mallows_loglik(\
        alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        prev_aug_rankings.slice(ii), metric\
      );
      double log_pcp = std::log(particle_correction_prob(ii));
      log_inc_wgt(ii) = loglik_1 - loglik_2 - log_pcp;
    }

    /* normalise weights ------------------------------------ */
    double maxw = arma::max(log_inc_wgt);
    arma::vec w = arma::exp(log_inc_wgt - maxw);
    arma::vec norm_wgt = w / arma::sum(w);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */
    /* Resample particles using multinomial resampling ------ */
    arma::uvec index = permutate_with_weights(norm_wgt, N);
    rho_samples.slice(tt + 1) = rho_samples.slice(tt + 1).rows(index);
    const arma::vec& asc = alpha_samples.col(tt + 1);
    alpha_samples.col(tt + 1) = asc.elem(index);
    aug_rankings = aug_rankings.slices(index);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (arma::uword ii = 0; ii < N; ++ii) {
      rho_samples.slice(tt + 1).row(ii) = metropolis_hastings_rho(\
        alpha_samples(ii, tt + 1), n_items, aug_rankings.slice(ii), metric,\
        rho_samples.slice(tt + 1).row(ii).t(), leap_size\
      ).t();
      alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(\
        alpha_samples(ii, tt + 1), n_items, aug_rankings.slice(ii), metric,\
          rho_samples.slice(tt + 1).row(ii).t(), logz_estimate,\
          alpha_prop_sd, lambda, alpha_max\
      );
      for (arma::uword jj = 0; jj < num_ranks; ++jj) {
        arma::vec mh_aug_result;
        if (aug_method == "random") {
          mh_aug_result = metropolis_hastings_aug_ranking(\
            alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(),\
            n_items, R_obs.slice(tt + 1).row(jj).t(),\
            aug_rankings.slice(ii).row(jj).t(), metric\
          );
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          mh_aug_result = metropolis_hastings_aug_ranking_pseudo(\
            alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(),\
            n_items, R_obs.slice(tt + 1).row(jj).t(),\
            aug_rankings.slice(ii).row(jj).t(), metric\
          );
        }
        aug_rankings.slice(ii).row(jj) = mh_aug_result.t();
      }
    }
  }

  /* ====================================================== */
  /* Post Processing                                        */
  /* ====================================================== */
  return Rcpp::List::create(
    Rcpp::Named("rho_samples") = rho_samples,
    Rcpp::Named("alpha_samples") = alpha_samples
  );
}
