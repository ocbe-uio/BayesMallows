#include "RcppArmadillo.h"
#include "distfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec get_is_estimate(arma::vec alpha_vector, int n,
                          std::string proposal_distribution,
                          std::string target_distribution,
                          int nmc
                          ) {

  // The dispersion parameter alpha comes as a vector value
  int num_alphas = alpha_vector.n_elem;

  arma::vec myind, ranks, prob, support, cpd, rho = arma::linspace<arma::vec>(1,n,n), u;

  arma::uvec inds;
  arma::vec logZ = arma::zeros(num_alphas);

  // Loop over the values of alpha
  for(int t = 0; t < num_alphas; ++t){
    // The current value of alpha
    double alpha = alpha_vector(t);
    // The current value of the partition function
    double Z = 0;

    // Loop over the Monte Carlo samples
    for(int i = 0; i < nmc; ++i){
      // Support set of the proposal distribution
      support = arma::linspace<arma::vec>(1,n,n);

      // Vector which holds the proposed ranks
      ranks = arma::zeros(n);

      // Probability of the ranks we get
      double q = 1;

      // n random uniform numbers
      u = arma::randu(n,1);

      // Loop over possible values given to item j in random order
      myind = arma::shuffle(arma::linspace(0, n-1, n));

      for(int j = 0; j < n; ++j){
        int jj = myind(j);
        prob = arma::zeros(n);
        // Find the elements that have not been taken yet
        inds = arma::find(support != 0);
        // Number of elements
        int k_max = inds.n_elem;

        // Reference vector
        arma::vec r1 = inds + arma::ones(k_max);
        // Sampled vector
        arma::vec r2 = rho(jj) * arma::ones(k_max);
        // Probability of sample
        prob(inds) = arma::exp(- alpha / n * arma::pow(arma::abs(r1 - r2), (proposal_distribution == "footrule") ? 1. : 2.));
        prob = prob / arma::accu(prob);
        cpd = arma::cumsum(prob);

        // Draw a random sampl
        inds = find(cpd > u(jj));
        ranks(jj) = inds(0) + 1;
        q = q * prob(inds(0));
        // Remove the realized rank from support
        support(ranks(jj) - 1) = 0;
      }
      // Increment the partition function
      Z += exp(- alpha / n * get_rank_distance(ranks, rho, target_distribution))/q;
    }
    // Average over the Monte Carlo samples
    logZ(t) = log(Z / nmc);
  }
  return logZ;
}
