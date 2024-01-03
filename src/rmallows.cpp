#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Sample from the Mallows distribution.
//'
//' Sample from the Mallows distribution with arbitrary distance metric using
//' a Metropolis-Hastings algorithm.
//'
//' @param rho0 Vector specifying the latent consensus ranking.
//' @param alpha0 Scalar specifying the scale parameter.
//' @param n_samples Integer specifying the number of random samples to generate.
//' @param burnin Integer specifying the number of iterations to discard as burn-in.
//' @param thinning Integer specifying the number of MCMC iterations to perform
//' between each time a random rank vector is sampled.
//' @param leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
//' @param metric Character string specifying the distance measure to use. Available
//' options are \code{"footrule"} (default), \code{"spearman"}, \code{"cayley"}, \code{"hamming"},
//' \code{"kendall"}, and \code{"ulam"}.
//'
//' @keywords internal
//'
//' @references \insertAllCited{}
//'
// [[Rcpp::export]]
arma::mat rmallows(
    arma::vec rho0,
    double alpha0,
    int n_samples,
    int burnin,
    int thinning,
    int leap_size = 1,
    std::string metric = "footrule"
  ) {
  auto distfun = choose_distance_function(metric);

  // The number of items ranked
  int n_items = rho0.n_elem;

  // Declare the matrix to hold the sampled ranks
  mat rho(n_items, n_samples);

  // Other variables used
  int rho_index = 0;

  // Vector to hold the iteration value of rho
  // Initializing it to the modal ranking
  vec rho_iter = rho0;

  int t = 1;
  // This is the Metropolis-Hastings loop
  while(rho_index < n_samples){
    // Check if the user has tried to stop the running
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    vec rho_proposal;
    uvec indices;
    double prob_backward, prob_forward;

    // Sample a proposal
    leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                   rho_iter, leap_size, distfun);

    distfun->update_leap_and_shift_indices(indices, n_items);

    // Compute the distances to current and proposed ranks
    const vec& rho1 = rho0(indices);
    double dist_new = distfun->d(rho1, rho_proposal(indices));
    double dist_old = distfun->d(rho1, rho_iter(indices));

    // Metropolis-Hastings ratio
    double ratio = - alpha0 / n_items * (dist_new - dist_old) +
      std::log(prob_backward) - std::log(prob_forward);

    // Draw a uniform random number
    double u = std::log(R::runif(0, 1));

    if(ratio > u){
      rho_iter = rho_proposal;
    }

    // Save every thinning'th iteration after burnin
    if((t > burnin) & (t % thinning == 0)){
      rho.col(rho_index) = rho_iter;
      ++rho_index;
    }

    ++t;
  }
  return rho;
}
