#include "RcppArmadillo.h"
#include "smc.h"
#include "partitionfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title SMC-Mallows new item rank (alpha fixed)
//' @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user
//' at each time step. Each correction and augmentation is done by filling in the missing item ranks randomly.
//' @param alpha A numeric value of the true scale parameter
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
//' @param alpha_prop_sd Numeric value specifying the standard deviation of the
//' lognormal proposal distribution used for \eqn{\alpha} in the
//' Metropolis-Hastings algorithm. Defaults to \code{0.1}.
//' @param lambda Strictly positive numeric value specifying the rate parameter
//' of the truncated exponential prior distribution of \eqn{\alpha}. Defaults
//' to \code{0.1}. When \code{n_cluster > 1}, each mixture component
//' \eqn{\alpha_{c}} has the same prior distribution.
//' @param alpha_max Maximum value of \code{alpha} in the truncated exponential
//' prior distribution.
//' @param aug_method A character string specifying the approach for filling in
//' the missing data, options are "pseudolikelihood" or "random".
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//' @return a 3d matrix containing the samples of rho and alpha from the SMC algorithm
//' @export
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_item_rank_alpha_fixed(
  double alpha,
  unsigned int& n_items,
  arma::cube& R_obs,
  std::string& metric,
  int& leap_size,
  unsigned int& N,
  unsigned int Time,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  int& mcmc_kernel_app,
  double alpha_prop_sd,
  double lambda,
  double alpha_max,
  std::string& aug_method,
  bool verbose = false
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */

  // Generate N initial samples of rho using the uniform prior
  arma::cube rho_samples(N, n_items, Time, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) {
    arma::uvec items_sample = arma::randperm(n_items) + 1;
    for (arma::uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  /* ====================================================== */
  /* Augment Rankings                                       */
  /* ====================================================== */
  unsigned int num_ranks = R_obs.n_rows;

  // each particle has its own set of augmented rankings
  arma::cube aug_rankings(num_ranks, n_items, N, arma::fill::zeros);
  arma::cube prev_aug_rankings(num_ranks, n_items, N, arma::fill::zeros);

  // augment incomplete ranks to initialise
  arma::ivec ranks = Rcpp::seq(1, n_items);

  // total correction prob
  arma::vec total_correction_prob = Rcpp::rep(1.0, N);

  // iterate through each observed ranking and create new "corrected" augmented rankings
  for (arma::uword ii = 0; ii < N; ++ii) {
    // set t-1 generation to old as we sample for t new
    prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

    for (arma::uword jj = 0; jj < num_ranks; ++jj) {
      arma::vec R_obs_slice_0_row_jj = R_obs.slice(0).row(jj).t();
      if (aug_method == "random") {
        // find elements missing from original observed ranking
        arma::vec partial_ranking = R_obs_slice_0_row_jj;

        Rcpp::NumericVector ranks_Cpp, partial_ranking_Cpp;
        ranks_Cpp = ranks;
        partial_ranking_Cpp = partial_ranking;
        Rcpp::NumericVector remaining_set = Rcpp::setdiff(ranks_Cpp, partial_ranking_Cpp);


        // create new agumented ranking by sampling remaining ranks from set uniformly
        arma::vec rset;
        int remaining_set_length = remaining_set.length();
        if (remaining_set_length == 1) {
          rset = Rcpp::as<arma::vec>(remaining_set);
        } else {
          rset = Rcpp::as<arma::vec>(Rcpp::sample(remaining_set, remaining_set.length()));
        }
        partial_ranking.elem(arma::find_nonfinite(partial_ranking)) = rset;


        aug_rankings.slice(ii).row(jj) = partial_ranking.t();
        // fill in missing ranks based on choice of augmentation method
        Rcpp::NumericVector remaining_set_length_Rcpp, remaining_set_length_Rcpp_fact;
        remaining_set_length_Rcpp = remaining_set_length;
        remaining_set_length_Rcpp_fact = Rcpp::factorial(remaining_set_length_Rcpp);
        double remaining_set_length_Rcpp_fact_dbl = Rcpp::as<double>(remaining_set_length_Rcpp_fact);
        total_correction_prob(ii) = total_correction_prob(ii) * (1 / remaining_set_length_Rcpp_fact_dbl);
      } else if ((aug_method == "pseudolikelihood") & ((metric == "footrule") | (metric == "spearman"))) {
        // find items missing from original observed ranking
        arma::uvec unranked_items = arma::find_nonfinite(R_obs_slice_0_row_jj);

        // find unallocated ranks from original observed ranking
        Rcpp::NumericVector rank_Rcpp, R_obs_slice_0_row_jj_Rcpp;
        R_obs_slice_0_row_jj_Rcpp = R_obs_slice_0_row_jj;
        rank_Rcpp = ranks;
        Rcpp::NumericVector remaining_set = Rcpp::setdiff(rank_Rcpp, R_obs_slice_0_row_jj_Rcpp);

        // randomly permute the unranked items to give the order in which they will be allocated
        arma::uvec item_ordering;
        item_ordering = arma::conv_to<arma::uvec>::from(arma::shuffle(unranked_items));
        Rcpp::List proposal = calculate_forward_probability(\
          item_ordering, R_obs_slice_0_row_jj, remaining_set, rho_samples.slice(0).row(ii).t(),\
          alpha, n_items, metric\
        );
        arma::vec a_rank = proposal["aug_ranking"];
        double f_prob = proposal["forward_prob"];
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
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue;

    /* Calculating log_z_alpha and log_likelihood ----------- */
    double log_z_alpha = get_partition_function(\
      n_items, alpha, cardinalities, logz_estimate, metric\
    );
    double log_likelihood = get_mallows_loglik(\
      alpha, rho_samples.slice(0).row(ii).t(), n_items,\
      aug_rankings.slice(ii), metric\
    ); // TODO: replace with log_lik_db? (#91)
    double log_tcp = std::log(total_correction_prob(ii));
    log_inc_wgt(ii) = log_likelihood - num_ranks * log_z_alpha - log_tcp;
  }

  /* normalise weights ------------------------------------ */
  double maxw = arma::max(log_inc_wgt);
  arma::vec w = arma::exp(log_inc_wgt - maxw);
  arma::vec norm_wgt = w / arma::sum(w);

  /* ====================================================== */
  /* Resample                                               */
  /* ====================================================== */
  /* Resample particles using multinomial resampling ------ */
  // Using norm_wgt_rcpp so that Rcpp::sample compiles. More details on
  // https://github.com/ocbe-uio/BayesMallows/issues/90#issuecomment-866614296
  Rcpp::NumericVector norm_wgt_rcpp;
  norm_wgt_rcpp = norm_wgt;
  arma::uvec index;
  index = Rcpp::as<arma::uvec>(Rcpp::sample(N, N, true, norm_wgt_rcpp));
  index -= 1;
  // index <- sample(1:N, prob = norm_wgt, size = N, replace = T)
  rho_samples.slice(0) = rho_samples.slice(0).rows(index);
  aug_rankings = aug_rankings.slices(index);

  /* ====================================================== */
  /* Move step                                              */
  /* ====================================================== */
  for (arma::uword ii = 0; ii < N; ++ii) {
    rho_samples.slice(0).row(ii) = metropolis_hastings_rho(\
      alpha, n_items, aug_rankings.slice(ii), metric,\
      rho_samples.slice(0).row(ii).t(), leap_size\
    ).t();
    for (arma::uword jj = 0; jj < num_ranks; ++jj) {
      arma::vec mh_aug_result;
      if (aug_method == "random") {
        mh_aug_result = metropolis_hastings_aug_ranking(\
          alpha, rho_samples.slice(0).row(ii).t(), n_items,\
          R_obs.slice(0).row(jj).t(), aug_rankings.slice(ii).row(jj).t(), metric\
        );
      } else if ((aug_method == "pseudolikelihood") & ((metric == "footrule") | (metric == "spearman"))) {
        mh_aug_result = metropolis_hastings_aug_ranking_pseudo(\
          alpha, rho_samples.slice(0).row(ii).t(), n_items,\
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
    if (verbose) REprintf("We are not on iteration %i out of %i \n", tt + 1, Time - 1);

    /* New Information -------------------------------------- */
    // new observed item ranks from each user, need to update augmented rankings
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);

    // total correction prob
    arma::vec total_correction_prob = Rcpp::rep(1.0, N);

    // iterate through each observed ranking and create new "corrected"
    // augmented rankings
    for (arma::uword ii = 0; ii < N; ++ii) {
      // set t-1 generation to old as we sample for t new
      prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

      // make the correction
      for (arma::uword jj = 0; jj < num_ranks; ++jj) {
        if (aug_method == "random") {
          Rcpp::List check_correction = correction_kernel(\
            R_obs.slice(tt + 1).row(jj).t(), aug_rankings.slice(ii).row(jj).t(),\
            n_items
          );
          arma::vec c_rank = check_correction["ranking"];
          double c_prob = check_correction["correction_prob"];
          aug_rankings.slice(ii).row(jj) = c_rank.t();
          total_correction_prob(ii) *= c_prob;
        } else if ((aug_method == "pseudolikelihood") & ((metric == "footrule") | (metric == "spearman"))) {
          Rcpp::List check_correction = correction_kernel_pseudo(\
            aug_rankings.slice(ii).row(jj).t(), R_obs.slice(tt + 1).row(jj).t(),\
            rho_samples.slice(tt + 1).row(ii).t(), alpha, n_items, metric\
          );
          arma::vec c_rank = check_correction["ranking"];
          double c_prob = check_correction["correction_prob"];
          aug_rankings.slice(ii).row(jj) = c_rank.t();
          // # these probs are in real scale
          total_correction_prob(ii) *= c_prob;
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
        alpha, rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        aug_rankings.slice(ii), metric\
      );
      double loglik_2 = get_mallows_loglik(\
        alpha, rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        prev_aug_rankings.slice(ii), metric\
      );
      double log_pcp = std::log(total_correction_prob(ii));
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
    // Using norm_wgt_rcpp so that Rcpp::sample compiles. More details on
    // https://github.com/ocbe-uio/BayesMallows/issues/90#issuecomment-866614296
    Rcpp::NumericVector norm_wgt_rcpp;
    norm_wgt_rcpp = norm_wgt;
    arma::uvec index;
    index = Rcpp::as<arma::uvec>(Rcpp::sample(N, N, true, norm_wgt_rcpp));
    index -= 1;
    rho_samples.slice(tt + 1) = rho_samples.slice(tt + 1).rows(index);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (arma::uword ii = 0; ii < N; ++ii) {
      rho_samples.slice(tt + 1).row(ii) = metropolis_hastings_rho(\
        alpha, n_items, aug_rankings.slice(ii), metric,\
        rho_samples.slice(tt + 1).row(ii).t(), leap_size\
      ).t();
      for (arma::uword jj = 0; jj < num_ranks; ++jj) {
        arma::vec mh_aug_result;
        if (aug_method == "random") {
          mh_aug_result = metropolis_hastings_aug_ranking(\
            alpha, rho_samples.slice(tt + 1).row(ii).t(),\
            n_items, R_obs.slice(tt + 1).row(jj).t(),\
            aug_rankings.slice(ii).row(jj).t(), metric\
          );
        } else if ((aug_method == "pseudolikelihood") & ((metric == "footrule") | (metric == "spearman"))) {
          mh_aug_result = metropolis_hastings_aug_ranking_pseudo(\
            alpha, rho_samples.slice(tt + 1).row(ii).t(),\
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
  return Rcpp::List::create((Rcpp::Named("rho_samples") = rho_samples));
}
