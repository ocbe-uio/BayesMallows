#include "RcppArmadillo.h"
#include "distances.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//' Compute importance sampling estimates of log partition function
//' for footrule and Spearman distances.
//'
//' @param alpha_vector Vector of alpha values at which to compute partition function.
//' @param n_items Integer specifying the number of ranked items.
//' @param metric Distance measure of the target Mallows distribution. Defaults to \code{footrule}.
//' @param nmc Number of Monte Carlo samples. Defaults to \code{1e4}.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec compute_importance_sampling_estimate(arma::vec alpha_vector, int n_items,
                          std::string metric = "footrule", int nmc = 1e4
                          ) {

  // The dispersion parameter alpha comes as a vector value
  int n_alphas = alpha_vector.n_elem;

  // The reference ranking
  arma::vec rho = arma::regspace<arma::vec>(1, n_items);

  // Vector which holds the result for all alphas
  arma::vec logZ = arma::zeros(n_alphas);

  // Loop over the values of alpha
  for(int t = 0; t < n_alphas; ++t){
    // The current value of alpha
    double alpha = alpha_vector(t);
    // The current value of the partition function
    arma::vec partfun(nmc);

    // Loop over the Monte Carlo samples
    for(int i = 0; i < nmc; ++i){
      // Support set of the proposal distribution
      arma::vec support = arma::regspace<arma::vec>(1, n_items);

      // Vector which holds the proposed ranks
      arma::vec ranks = arma::zeros(n_items);
      arma::vec ranks2 = arma::zeros(n_items);

      // Probability of the ranks we get
      double log_q = 0;

      // n_items random uniform numbers
      arma::vec u = arma::log(arma::randu(n_items));

      // Loop over possible values given to item j in random order
      arma::vec myind = arma::shuffle(arma::regspace(0, n_items - 1));

      for(int j = 0; j < n_items; ++j){
        int jj = myind(j);
        // Find the elements that have not been taken yet
        arma::uvec inds = arma::find(support != 0);
        arma::vec log_prob(inds.size());

        // Number of elements
        int k_max = inds.n_elem;

        // Reference vector
        arma::vec r1 = inds + arma::ones(k_max);
        // Sampled vector
        arma::vec r2 = rho(jj) * arma::ones(k_max);
        // Probability of sample. Note that this is a vector quantity.
        log_prob = - alpha / n_items * arma::pow(arma::abs(r1 - r2), (metric == "footrule") ? 1. : 2.);
        log_prob = log_prob - std::log(arma::accu(arma::exp(log_prob)));
        arma::vec log_cpd = arma::log(arma::cumsum(arma::exp(log_prob)));

        // Draw a random sample
        int item_index = arma::as_scalar(arma::find(log_cpd > u(jj), 1));
        ranks(jj) = arma::as_scalar(inds(item_index)) + 1;

        log_q += log_prob(item_index);

        // Remove the realized rank from support
        support(ranks(jj) - 1) = 0;
      }

      // Increment the partition function
      partfun(i) = - alpha / n_items * get_rank_distance(ranks, rho, metric) - log_q;
    }
    // Average over the Monte Carlo samples
    // Using this trick: https://www.xarg.org/2016/06/the-log-sum-exp-trick-in-machine-learning/
    double maxval = arma::max(partfun);
    logZ(t) = maxval + std::log(arma::accu(arma::exp(partfun - maxval))) - std::log(static_cast<double>(nmc));
  }
  return logZ;
}
