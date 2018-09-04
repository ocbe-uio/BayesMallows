#include "RcppArmadillo.h"
#include "leapandshift.h"
#include "distfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Worker function for sampling from the Mallows distribution.
//'
//' @param rho0 Vector specifying the latent consensus ranking.
//' @param alpha0 Scalar specifying the scale parameter.
//' @param n_samples Integer specyfing the number of random samples to generate.
//' @param burnin Integer specifying the number of iterations to discard as burn-in.
//' @param thinning Integer specifying the number of MCMC iterations to perform
//' between each time a random rank vector is sampled.
//' @param leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
//' @param metric Character string specifying the distance measure to use. Available
//' options are \code{"footrule"} (default) and \code{"Spearman"}.
//'
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec rmallows(
    arma::vec rho0,
    double alpha0,
    int n_samples,
    int burnin,
    int thinning,
    int leap_size = 1,
    std::string metric = "footrule"
  ){

  // The number of items ranked
  int n_items = rho0.n_elem;

  // Declare the matrix to hold the sampled ranks
  arma::mat rho(n_items, n_samples);

  // Other variables used
  int rho_index = 0;

  // Vector to hold the iteration value of rho
  // Initializing it to the modal ranking
  arma::vec rho_iter = rho0;

  int t = 1;
  // This is the Metropolis-Hastings loop
  while((t < burnin) | (rho_index < n_samples)){

    // Check if the user has tried to stop the running
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    arma::vec rho_proposal;
    arma::uvec indices;
    double prob_backward, prob_forward;

    // Sample a proposal
    leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                   rho_iter, leap_size);

    // Metropolis-Hastings ratio
    double ratio = - alpha0 / n_items * get_rank_distance(rho_proposal, rho0, metric) +
      log(prob_backward) + alpha0 / n_items * get_rank_distance(rho_iter, rho0, metric) -
      log(prob_forward);

    // Draw a uniform random number
    double u = log(arma::randu<double>());

    if(ratio > u){
      rho_iter = rho_proposal;
    }

    // Save every thinning'th iteration after burnin
    if(t > burnin & t % thinning == 0){
      rho.col(rho_index) = rho_iter;
      ++rho_index;
    }

    ++t;
  }

  return rho;


}



