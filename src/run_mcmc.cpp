#include <math.h>
#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// We put all the function declarations here for clarity
int binomial_coefficient(int n, int k);
double get_rank_distance(arma::vec, arma::vec, std::string);
double get_partition_function(int, double, arma::vec, std::string);
arma::vec get_summation_distances(int, arma::vec, std::string);
Rcpp::List leap_and_shift(arma::vec, int);
double rank_dist_matrix(arma::mat, arma::vec, std::string);
int factorial(int);

//' Worker function for computing the posterior distribtuion.
//'
//' @param R A set of complete rankings, with one sample per column.
//' With N samples and n items, R is n x N.
//' @param nmc Number of Monte Carlo samples.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function. Otherwise,
//' please provided an arbitrary vector.
//' @param metric The distance metric to use. On of \code{"spearman"},
//' \code{"footrule"}, \code{"kendall"}, \code{"cayley"}, or
//' \code{"hamming"}.
//' @param L Leap-and-shift step size.
//' @param sd_alpha Standard deviation of proposal distribution for alpha.
//' @param alpha_init Initial value of alpha.
//' @param alpha_jump How many times should we sample \code{rho} between
//' each time we sample \code{alpha}. Setting \code{alpha_jump} to a high
//' number can significantly speed up computation time, since we then do not
//' have to do expensive computation of the partition function.
//' @param lambda Parameter of the prior distribution.
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat R, int nmc,
                    arma::vec cardinalities,
                    std::string metric = "footrule",
                    int L = 1, double sd_alpha = 0.5,
                    double alpha_init = 5, int alpha_jump = 1,
                    double lambda = 0.1){

  // The number of items ranked
  int n = R.n_rows;

  // The number of assessors
  int N = R.n_cols;

  // Number of alpha values to store.
  int n_alpha = ceil(nmc * 1.0 / alpha_jump);

  // Declare the matrix to hold the latent ranks
  // Note: Armadillo matrices are stored in column-major ordering. Hence,
  // we put the items along the column, since they are going to be accessed at the
  // same time for a given Monte Carlo sample.
  arma::mat rho(n, nmc);

  // Set the initial latent rank value
  rho.col(0) = arma::linspace<arma::vec>(1, n, n);

  // Declare the vector to hold the scaling parameter alpha
  arma::vec alpha(n_alpha);

  // Set the initial alpha value
  alpha(0) = alpha_init;

  // Declare indicator vectors to hold acceptance or not
  arma::vec alpha_acceptance(n_alpha), rho_acceptance(nmc);

  // Set the initial values;
  alpha_acceptance(0) = 1;
  rho_acceptance(0) = 1;

  // Metropolis-Hastings ratio
  double ratio;

  // Uniform random number which we will draw in Metropolis-Hastings algorithm
  double u;

  // Other variables used
  int alpha_index = 0;
  double alpha_old = alpha(0);

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0
  for(int t = 1; t < nmc; ++t){

    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    // Save current parameter values
    arma::vec rho_old = rho.col(t - 1);

    if(t % alpha_jump == 0) {

      // Increment to the index we are going to update
      ++alpha_index;

      // Sample an alpha proposal (normal on the log scale)
      double alpha_proposal = exp(arma::randn<double>() * sd_alpha +
                                  log(alpha(alpha_index - 1)));

      double rank_dist = rank_dist_matrix(R, rho_old, metric);

      // Difference between current and proposed alpha
      double alpha_diff = alpha(alpha_index - 1) - alpha_proposal;

      // Compute the Metropolis-Hastings ratio
      ratio = (alpha(alpha_index - 1) - alpha_proposal) / n * rank_dist +
        lambda * alpha_diff +
        N * (
            get_partition_function(n, alpha(alpha_index - 1), cardinalities, metric) -
              get_partition_function(n, alpha_proposal, cardinalities, metric)
        ) + log(alpha_proposal) - log(alpha(alpha_index - 1));

      // Draw a uniform random number
      u = log(arma::randu<double>());


      if(ratio > u){
        alpha(alpha_index) = alpha_proposal;
        alpha_acceptance(alpha_index) = 1;
      } else {
        alpha(alpha_index) = alpha(alpha_index - 1);
        alpha_acceptance(alpha_index) = 0;
      }

      // Next time around, this is alpha_old
      // The sampling of rho below should also
      // use this alpha_old, which then in the first
      // step is actually the newest alpha
      alpha_old = alpha(alpha_index);
    }

    // Sample a rank proposal
    Rcpp::List ls_proposal = leap_and_shift(rho_old, L);

    // Save some of the variables
    arma::vec rho_proposal = ls_proposal["proposal"];
    arma::uvec indices = ls_proposal["indices"];
    double prob_backward = ls_proposal["prob_backward"];
    double prob_forward = ls_proposal["prob_forward"];

    // Compute the distances to current and proposed ranks
    double dist_new = rank_dist_matrix(R.rows(indices), rho_proposal(indices), metric);
    double dist_old = rank_dist_matrix(R.rows(indices), rho_old(indices), metric);

    // Metropolis-Hastings ratio
    ratio = - alpha_old / n * (dist_new - dist_old) +
      log(prob_backward) - log(prob_forward);

    // Draw a uniform random number
    u = log(arma::randu<double>());

    if(ratio > u){
      rho.col(t) = rho_proposal;
      rho_acceptance(t) = 1;
    } else {
      rho.col(t) = rho.col(t - 1);
      rho_acceptance(t) = 0;
    }


  }

  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance
  );
}

