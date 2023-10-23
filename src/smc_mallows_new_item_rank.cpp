#include <RcppArmadillo.h>
#include "smc.h"
#include "sample.h"
#include "parameterupdates.h"
#include "setdiff.h"
#include "misc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void new_items_move_step(
    cube& rho_samples,
    mat& alpha_samples,
    cube& aug_rankings,
    const cube& rankings,
    const std::string& aug_method,
    const Rcpp::Nullable<arma::vec>& logz_estimate,
    const Rcpp::Nullable<arma::vec>& cardinalities,
    const double& alpha,
    const uword& ttplus1,
    const bool& alpha_fixed,
    const std::string& metric = "footrule",
    const int& leap_size = 1,
    const double& alpha_prop_sd = 0.5,
    const double& alpha_max = 1e6,
    const double& lambda = 0.1
){
  const uword num_ranks = rankings.n_rows;
  const uword n_particles = rho_samples.n_rows;
  int n_items = rho_samples.n_cols;
  for (uword ii = 0; ii < n_particles; ++ii) {
    rho_samples.slice(ttplus1).row(ii) = metropolis_hastings_rho(
      alpha_fixed ? alpha : alpha_samples(ii, ttplus1), n_items, aug_rankings.slice(ii),
      rho_samples.slice(ttplus1).row(ii).t(), metric, leap_size
    ).t();
    if(!alpha_fixed){
      alpha_samples(ii, ttplus1) = metropolis_hastings_alpha(
        alpha_samples(ii, ttplus1), n_items, aug_rankings.slice(ii),
        rho_samples.slice(ttplus1).row(ii).t(), logz_estimate,
        cardinalities,
        metric, alpha_prop_sd, alpha_max, lambda
      );
    }

    for (uword jj = 0; jj < num_ranks; ++jj) {
      vec mh_aug_result = metropolis_hastings_aug_ranking(                                  \
          alpha_fixed ? alpha : alpha_samples(ii, ttplus1), rho_samples.slice(ttplus1).row(ii).t(),\
          n_items, rankings.slice(ttplus1).row(jj).t(),                                              \
          aug_rankings.slice(ii).row(jj).t(),
          is_pseudo(aug_method, metric),
          metric
        );
      aug_rankings.slice(ii).row(jj) = mh_aug_result.t();
    }
  }
}

arma::cube augment_rankings(
  const unsigned int& n_items,
  arma::cube& rankings,
  const unsigned int& n_particles,
  arma::cube rho_samples,
  arma::mat alpha_samples,
  const double alpha = 0,
  const std::string& aug_method = "random",
  const bool alpha_fixed = false,
  const std::string& metric = "footrule"
) {
  /* ====================================================== */
  /* Augment Rankings                                       */
  /* ====================================================== */
  const unsigned int& num_ranks = rankings.n_rows;

  // each particle has its own set of augmented rankings
  cube aug_rankings(num_ranks, n_items, n_particles, fill::zeros);
  cube prev_aug_rankings(num_ranks, n_items, n_particles, fill::zeros);

  // augment incomplete ranks to initialise
  const ivec ranks = regspace<ivec>(1, n_items);

  // total correction prob
  vec total_correction_prob = ones(n_particles);

  // iterate through each observed ranking and create new "corrected" augmented
  // rankings
  for (uword ii = 0; ii < n_particles; ++ii) {
    // set t-1 generation to old as we sample for t new
    prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

    // make the correction
    for (uword jj = 0; jj < num_ranks; ++jj) {
      // fill in missing ranks based on choice of augmentation method
      vec rankings_slice_0_row_jj = rankings.slice(0).row(jj).t();
      const vec remaining_set = setdiff_template(ranks, rankings_slice_0_row_jj);
      const bool pseudo = is_pseudo(aug_method, metric);
      if (!pseudo) {
        // create new augmented ranking by sampling remaining ranks from set uniformly
        vec rset = shuffle(remaining_set);

        vec partial_ranking = rankings_slice_0_row_jj;
        partial_ranking.elem(find_nonfinite(partial_ranking)) = rset;

        aug_rankings.slice(ii).row(jj) = partial_ranking.t();
      } else {
        // find items missing from original observed ranking
        const uvec& unranked_items = find_nonfinite(rankings_slice_0_row_jj);
        // randomly permute the unranked items to give the order in which they will be allocated
        uvec item_ordering = shuffle(unranked_items);
        const Rcpp::List proposal = calculate_forward_probability(\
          item_ordering, rankings_slice_0_row_jj, remaining_set, rho_samples.slice(0).row(ii).t(),\
          alpha_fixed ? alpha : alpha_samples(ii, 0), n_items, metric\
        );
        const vec& a_rank = proposal["aug_ranking"];
        aug_rankings(span(jj), span::all, span(ii)) = a_rank;
      }
    }
  }
  return(aug_rankings);
}

//' @title SMC-Mallows new item rank
//' @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user
//' at each time step. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.
//' @param n_items Integer is the number of items in a ranking
//' @param rankings 3D matrix of size n_assessors by n_items by timesteps containing a set of observed rankings of timesteps time steps
//' @param metric A character string specifying the distance metric to use in the
//' Bayesian Mallows Model. Available options are \code{"footrule"},
//' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//' \code{"ulam"}.
//' @param leap_size leap_size Integer specifying the step size of the leap-and-shift
//' proposal distribution
//' @param n_particles Integer specifying the number of particles
//' @param timesteps Integer specifying the number of time steps in the SMC algorithm
//' @param logz_estimate Estimate of the partition function, computed with
//' \code{\link{estimate_partition_function}}.
//' @param cardinalities Cardinalities for exact computation of partition function,
//' returned from \code{\link{prepare_partition_function}}.
//' @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel
//' @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha
//' @param lambda Strictly positive numeric value specifying the rate parameter
//' of the truncated exponential prior distribution of alpha.
//' @param alpha_max  Maximum value of alpha in the truncated exponential
//' prior distribution.
//' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
//' @param verbose Logical specifying whether to print out the progress of the
//' SMC-Mallows algorithm. Defaults to \code{FALSE}.
//' @param alpha_fixed Logical indicating whether to sample \code{alpha} or not.
//' @param alpha numeric value of the scale parameter.
//' @param aug_rankings_init Initial values for augmented rankings.
//' @param rho_samples_init Initial values for rho samples.
//' @param alpha_samples_init Initial values for alpha samples.
//'
//' @return a 3d matrix containing: the samples of: rho, alpha and the augmented rankings, and the effective sample size at each iteration of the SMC algorithm.
//'
//' @noRd
//'
//' @family modeling
//'
// [[Rcpp::export]]
Rcpp::List smc_mallows_new_item_rank_cpp(
  const unsigned int& n_items,
  arma::cube& rankings,
  const unsigned int& n_particles,
  const unsigned int timesteps,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const Rcpp::Nullable<arma::vec> cardinalities,
  const int& mcmc_kernel_app,
  Rcpp::Nullable<arma::cube> aug_rankings_init = R_NilValue,
  Rcpp::Nullable<arma::mat> rho_samples_init = R_NilValue,
  arma::vec alpha_samples_init = 0,
  const double alpha = 0,
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1,
  const double alpha_max = 1e6,
  const std::string& aug_method = "random",
  const bool verbose = false,
  const bool alpha_fixed = false,
  const std::string& metric = "footrule",
  const int& leap_size = 1
) {
  /* ====================================================== */
  /* Initialise Phase                                       */
  /* ====================================================== */

  // Generate n_particles initial samples of rho using the uniform prior
  cube rho_samples(n_particles, n_items, timesteps);
  if (rho_samples_init.isNotNull()) {
    rho_samples.slice(0) = Rcpp::as<arma::mat>(rho_samples_init);
  } else {
    rho_samples.slice(0) = initialize_rho(n_items, n_particles).t();
  }

  mat alpha_samples;
  if(!alpha_fixed){
    // If alpha_fixed = false, alpha_samples needs to be generated from
    // alpha_samples_init
    alpha_samples = zeros(n_particles, timesteps);
    if (alpha_samples_init.n_elem != n_particles) {
      alpha_samples_init = initialize_alpha(n_particles);
    }
    alpha_samples.col(0) = alpha_samples_init;
  }

  /* generate vector to store ESS */
  rowvec ESS_vec(timesteps);
  ESS_vec(0) = 1;

  /* ====================================================== */
  /* Augment Rankings                                       */
  /* ====================================================== */
  const unsigned int& num_ranks = rankings.n_rows;

  // each particle has its own set of augmented rankings
  cube aug_rankings(num_ranks, n_items, n_particles, fill::zeros);
  cube prev_aug_rankings(num_ranks, n_items, n_particles, fill::zeros);

  // augment incomplete ranks to initialise
  const ivec ranks = regspace<ivec>(1, n_items);

  // augment rankings for proposal
  arma::cube aug_rankings_init_2 = aug_rankings_init.isNotNull() ? Rcpp::as<arma::cube>(aug_rankings_init) : arma::cube(0,0,0);
  if (aug_rankings_init_2.size() == 0) {
    aug_rankings = augment_rankings(
      n_items, rankings, n_particles, rho_samples, alpha_samples, alpha, aug_method,
      alpha_fixed, metric
    );
  } else {
      aug_rankings = aug_rankings_init_2;
  }
  // set the old augmentation as we will compare the t^th augmentation
  // to the (t-1)^th augmentation
  prev_aug_rankings = aug_rankings;

  /* ====================================================== */
  /* Loop for t=1,...,timesteps                                  */
  /* ====================================================== */
  // Here, we attempt the SMC sampler
  const bool pseudo = is_pseudo(aug_method, metric);
  for (uword tt = 0; tt < timesteps - 1; ++tt) {
    if (verbose) REprintf("iteration %i out of %i \n", tt + 1, timesteps - 1);
    /* New Information -------------------------------------- */
    // new observed item ranks from each user, need to update augmented rankings
    rho_samples.slice(tt + 1) = rho_samples.slice(tt);
    if(!alpha_fixed) alpha_samples.col(tt + 1) = alpha_samples.col(tt);

    // total correction prob
    vec particle_correction_prob = ones(n_particles);

    // iterate through each observed ranking and create new "corrected"
    // augmented rankings
    for (uword ii = 0; ii < n_particles; ++ii) {
      // set t-1 generation to old as we sample for t new
      prev_aug_rankings.slice(ii) = aug_rankings.slice(ii);

      // make the correction
      for (uword jj = 0; jj < num_ranks; ++jj) {
        Rcpp::List check_correction;
        if (!pseudo) {
          check_correction = correction_kernel(
            rankings.slice(tt + 1).row(jj).t(), aug_rankings.slice(ii).row(jj).t(),
            n_items
          );
        } else {
          check_correction = correction_kernel_pseudo(
            aug_rankings.slice(ii).row(jj).t(), rankings.slice(tt + 1).row(jj).t(),
            rho_samples.slice(tt + 1).row(ii).t(),
            alpha_fixed ? alpha : alpha_samples(ii, tt + 1), n_items, metric
          );
        }
        const vec& c_rank = check_correction["ranking"];
        aug_rankings(span(jj), span::all, span(ii)) = c_rank;
        const double& c_prob = check_correction["correction_prob"];
        particle_correction_prob(ii) *= c_prob;
      }
    }
    /* ====================================================== */
    /* Re-weight                                              */
    /* ====================================================== */

    // incremental weight for each particle, based on new observed rankings
    vec log_inc_wgt(n_particles, fill::zeros);

    for (uword ii = 0; ii < n_particles; ++ii) {
      // evaluate the log estimate of the partition function for a particular
      // value of alpha

      /* Calculating log_z_alpha and log_likelihood ----------- */
      const double& loglik_1 = get_exponent_sum(\
        alpha_fixed ? alpha : alpha_samples(ii, tt + 1),\
        rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        aug_rankings.slice(ii), metric\
      );
      const double& loglik_2 = get_exponent_sum(\
        alpha_fixed ? alpha : alpha_samples(ii, tt + 1),\
        rho_samples.slice(tt + 1).row(ii).t(), n_items,\
        prev_aug_rankings.slice(ii), metric\
      );
      const double& log_pcp = std::log(particle_correction_prob(ii));
      log_inc_wgt(ii) = loglik_1 - loglik_2 - log_pcp;
    }

    // update weights
    vec norm_wgt = normalize_weights(log_inc_wgt);

    /* store ESS = sum(w)^2/sum(w^2) */
    ESS_vec(tt + 1) = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);

    /* ====================================================== */
    /* Resample                                               */
    /* ====================================================== */
    /* Resample particles using multinomial resampling ------ */
    uvec index = sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);
    rho_samples.slice(tt + 1) = rho_samples.slice(tt + 1).rows(index);
    if(!alpha_fixed){
      const vec& asc = alpha_samples.col(tt + 1);
      alpha_samples.col(tt + 1) = asc.elem(index);
    }

    aug_rankings = aug_rankings.slices(index);

    /* ====================================================== */
    /* Move step                                              */
    /* ====================================================== */
    new_items_move_step(
      rho_samples, alpha_samples, aug_rankings, rankings, aug_method,
      logz_estimate, cardinalities, alpha, tt + 1,
      alpha_fixed, metric, leap_size, alpha_prop_sd, alpha_max, lambda);
  }

  /* ====================================================== */
  /* Post Processing                                        */
  /* ====================================================== */
  if (alpha_fixed) {
    Rcpp::List particle_history = Rcpp::List::create(
      Rcpp::Named("rho_samples") = rho_samples,
      Rcpp::Named("augmented_rankings") = aug_rankings,
      Rcpp::Named("ESS") = ESS_vec
    );
    particle_history.attr("class") = "SMCMallows";
    return particle_history;
  } else {
    Rcpp::List particle_history = Rcpp::List::create(
      Rcpp::Named("rho_samples") = rho_samples,
      Rcpp::Named("alpha_samples") = alpha_samples,
      Rcpp::Named("augmented_rankings") = aug_rankings,
      Rcpp::Named("ESS") = ESS_vec
    );
    particle_history.attr("class") = "SMCMallows";
    return particle_history;
  }
}
