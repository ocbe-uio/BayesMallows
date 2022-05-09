#include <RcppArmadillo.h>
#include "misc.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"

using namespace arma;

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
//' @return a set of particles each containing the values of rho and alpha and the effective sample size (ESS) at each iteration of the SMC
//' algorithm as well as the set of augmented rankings at the final iteration.
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
  cube rho_samples(N, n_items, Time, fill::zeros);
  for (uword i = 0; i < N; ++i) {
    const uvec items_sample = randperm(n_items) + 1;
    for (uword j = 0; j < n_items; ++j) {
      rho_samples(i, j, 0) = items_sample(j);
    }
  }

  /* generate alpha samples using exponential prior ------- */
  mat alpha_samples(N, Time);
  const vec alpha_samples_0 = Rcpp::rexp(N, 1);
  alpha_samples.col(0) = alpha_samples_0;

  /* ====================================================== */
  /* Augment Rankings                                       */
  /* ====================================================== */
  const unsigned int& num_ranks = R_obs.n_rows;

  // each particle has its own set of augmented rankings
  cube aug_rankings(num_ranks, n_items, N, fill::zeros);
  cube prev_aug_rankings(num_ranks, n_items, N, fill::zeros);

  // augment incomplete ranks to initialise
  const ivec ranks = Rcpp::seq(1, n_items);

  // total correction prob
  vec total_correction_prob = Rcpp::rep(1.0, N);

  // iterate through each observed ranking and create new "corrected" augmented rankings
  for (uword ii = 0; ii < N; ++ii) {
    // set t-1 generation to old as we sample for t new
    prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

    // make the correction
    for (uword jj = 0; jj < num_ranks; ++jj) {
      // fill in missing ranks based on choice of augmentation method
      vec R_obs_slice_0_row_jj = R_obs.slice(0).row(jj).t();
      const vec remaining_set = setdiff_template(ranks, R_obs_slice_0_row_jj);
      if (aug_method == "random") {
        // create new augmented ranking by sampling remaining ranks from set uniformly
        vec rset = shuffle(remaining_set);

        vec partial_ranking = R_obs_slice_0_row_jj;
        partial_ranking.elem(find_nonfinite(partial_ranking)) = rset;

        aug_rankings.slice(ii).row(jj) = partial_ranking.t();
        total_correction_prob(ii) = divide_by_fact(total_correction_prob(ii), remaining_set.size());
      } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
        // find items missing from original observed ranking
        const uvec& unranked_items = find_nonfinite(R_obs_slice_0_row_jj);
        // randomly permute the unranked items to give the order in which they will be allocated
        uvec item_ordering;
        item_ordering = conv_to<uvec>::from(shuffle(unranked_items));
        const Rcpp::List proposal = calculate_forward_probability(\
          item_ordering, R_obs_slice_0_row_jj, remaining_set, rho_samples.slice(0).row(ii).t(),\
          alpha_samples(ii, 0), n_items, metric\
        );
        const vec& a_rank = proposal["aug_ranking"];
        const double& f_prob = proposal["forward_prob"];
        aug_rankings(span(jj), span::all, span(ii)) = a_rank;
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
  vec log_inc_wgt(N, fill::zeros);

  for (uword ii = 0; ii < N; ++ii) {
    // evaluate the log estimate of the partition function for a particular
    // value of alpha

    /* Initializing variables ------------------------------- */
    const Rcpp::Nullable<vec>& cardinalities = R_NilValue;

    /* Calculating log_z_alpha and log_likelihood ----------- */
    const double& log_z_alpha = get_partition_function(\
      n_items, alpha_samples(ii, 0), cardinalities, logz_estimate, metric\
    );
    const double& log_likelihood = get_exponent_sum(\
      alpha_samples(ii, 0), rho_samples.slice(0).row(ii).t(), n_items,\
      aug_rankings.slice(ii), metric\
    );
    const double& log_tcp = std::log(total_correction_prob(ii));
    log_inc_wgt(ii) = log_likelihood - num_ranks * log_z_alpha - log_tcp;
  }

  /* normalise weights ------------------------------------ */
  const double& maxw = max(log_inc_wgt);
  const vec& w = exp(log_inc_wgt - maxw);
  const vec& norm_wgt = w / sum(w);

  /* store ESS result in first entry of the vector ESS_vec */
  /* generate vector to store ESS */
  double ess_init = sum(norm_wgt) * sum(norm_wgt) / sum(norm_wgt % norm_wgt);
  rowvec ESS_vec(Time, fill::zeros);
  ESS_vec(0) = ess_init;

  /* ====================================================== */
  /* Resample                                               */
  /* ====================================================== */
  /* Resample particles using multinomial resampling ------ */
  uvec index = permute_with_weights(norm_wgt, N);
  rho_samples.slice(0) = rho_samples.slice(0).rows(index);
  const vec& asc = alpha_samples.col(0);
  alpha_samples.col(0) = asc.elem(index);
  aug_rankings = aug_rankings.slices(index);

  /* ====================================================== */
  /* Move step                                              */
  /* ====================================================== */
  for (uword ii = 0; ii < N; ++ii) {
    rho_samples.slice(0).row(ii) = metropolis_hastings_rho(\
      alpha_samples(ii, 0), n_items, aug_rankings.slice(ii), metric,\
      rho_samples.slice(0).row(ii).t(), leap_size\
    ).t();
    alpha_samples(ii, 0) = metropolis_hastings_alpha(\
      alpha_samples(ii, 0), n_items, aug_rankings.slice(ii), metric,\
        rho_samples.slice(0).row(ii).t(), logz_estimate,\
        alpha_prop_sd, lambda, alpha_max\
    );
    for (uword jj = 0; jj < num_ranks; ++jj) {
      vec mh_aug_result;
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

  for (uword tt = 0; tt < Time - 1; ++tt) {
    if (verbose) REprintf("iteration %i out of %i \n", tt + 1, Time - 1);

    /* New Information -------------------------------------- */
    // new observed item ranks from each user, need to update augmented rankings
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    alpha_samples.col(tt + 1) = alpha_samples.col(tt);

    // total correction prob
    vec particle_correction_prob = Rcpp::rep(1.0, N);

    // iterate through each observed ranking and create new "corrected"
    // augmented rankings

    for (uword ii = 0; ii < N; ++ii) {
      // set t-1 generation to old as we sample for t new
      prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

      // make the correction
      for (uword jj = 0; jj < num_ranks; ++jj) {
        if (aug_method == "random") {
          const Rcpp::List check_correction = correction_kernel(\
            R_obs.slice(tt + 1).row(jj).t(), aug_rankings.slice(ii).row(jj).t(),\
            n_items
          );
          const vec c_rank = check_correction["ranking"];
          const double c_prob = check_correction["correction_prob"];
          aug_rankings.slice(ii).row(jj) = c_rank.t();
          particle_correction_prob(ii) *= c_prob;
        } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
          const Rcpp::List check_correction = correction_kernel_pseudo(\
            aug_rankings.slice(ii).row(jj).t(), R_obs.slice(tt + 1).row(jj).t(),\
            rho_samples.slice(tt + 1).row(ii).t(), alpha_samples(ii, tt + 1),\
            n_items, metric
          );
          const vec c_rank = check_correction["ranking"];
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
    vec log_inc_wgt(N, fill::zeros);
    for (uword ii = 0; ii < N; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Calculating log_z_alpha and log_likelihood ----------- */
      double loglik_1 = get_exponent_sum(\
        alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        aug_rankings.slice(ii), metric\
      );
      double loglik_2 = get_exponent_sum(\
        alpha_samples(ii, tt + 1), rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        prev_aug_rankings.slice(ii), metric\
      );
      double log_pcp = std::log(particle_correction_prob(ii));
      log_inc_wgt(ii) = loglik_1 - loglik_2 - log_pcp;
    }

    /* normalise weights ------------------------------------ */
    double maxw = max(log_inc_wgt);
    vec w = exp(log_inc_wgt - maxw);
    vec norm_wgt = w / sum(w);
    /* store ESS = sum(w)^2/sum(w^2) */
    ESS_vec(tt + 1) = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */
    /* Resample particles using multinomial resampling ------ */
    uvec index = permute_with_weights(norm_wgt, N);
    rho_samples.slice(tt + 1) = rho_samples.slice(tt + 1).rows(index);
    const vec& asc = alpha_samples.col(tt + 1);
    alpha_samples.col(tt + 1) = asc.elem(index);
    aug_rankings = aug_rankings.slices(index);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    for (uword ii = 0; ii < N; ++ii) {
      rho_samples.slice(tt + 1).row(ii) = metropolis_hastings_rho(\
        alpha_samples(ii, tt + 1), n_items, aug_rankings.slice(ii), metric,\
        rho_samples.slice(tt + 1).row(ii).t(), leap_size\
      ).t();
      alpha_samples(ii, tt + 1) = metropolis_hastings_alpha(\
        alpha_samples(ii, tt + 1), n_items, aug_rankings.slice(ii), metric,\
          rho_samples.slice(tt + 1).row(ii).t(), logz_estimate,\
          alpha_prop_sd, lambda, alpha_max\
      );
      for (uword jj = 0; jj < num_ranks; ++jj) {
        vec mh_aug_result;
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
    Rcpp::Named("alpha_samples") = alpha_samples,
    Rcpp::Named("augmented_rankings") = aug_rankings,
    Rcpp::Named("ESS") = ESS_vec
  );
}
