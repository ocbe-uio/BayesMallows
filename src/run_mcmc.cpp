#include <math.h>
#include "RcppArmadillo.h"
#include "parameterupdates.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Worker function for computing the posterior distribtuion.
//'
//' @param R A set of complete rankings, with one sample per column.
//' With N samples and n items, R is n x N.
//' @param nmc Number of Monte Carlo samples.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function. Defaults to
//' \code{R_NilValue}.
//' @param metric The distance metric to use. One of \code{"spearman"},
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
                    Rcpp::Nullable<arma::vec> cardinalities,
                    Rcpp::Nullable<arma::vec> is_fit,
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

  // Other variables used
  int alpha_index = 0;
  double alpha_old = alpha(0);
  arma::vec rho_old = rho.col(0);

  // This is the Metropolis-Hastings loop
  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0
  for(int t = 1; t < nmc; ++t){

    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    if(t % alpha_jump == 0) {
      // Call the void function which updates alpha by reference
      update_alpha(alpha, alpha_acceptance, alpha_old, R, alpha_index,
        rho_old, sd_alpha, metric, lambda, n, N, cardinalities, is_fit);
    }

    // Call the void function which updates rho by reference
    update_rho(rho, rho_acceptance, rho_old, alpha_old, L, R, metric, n, t);

  }

  // Finally return the MCMC samples
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance
  );
}



