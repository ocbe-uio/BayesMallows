#include <RcppArmadillo.h>
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec compute_importance_sampling_estimate(
    arma::vec alpha_vector, int n_items,
    std::string metric = "footrule", int nmc = 1e4
                          ) {
  auto distfun = choose_distance_function(metric);
  // The dispersion parameter alpha comes as a vector value
  int n_alphas = alpha_vector.n_elem;

  // The reference ranking
  vec rho = regspace<vec>(1, n_items);

  // Vector which holds the result for all alphas
  vec logZ = zeros(n_alphas);

  // Loop over the values of alpha
  for(int t = 0; t < n_alphas; ++t){
    // The current value of alpha
    double alpha = alpha_vector(t);
    // The current value of the partition function
    vec partfun(nmc);

    // Loop over the Monte Carlo samples
    for(int i = 0; i < nmc; ++i){
      // Support set of the proposal distribution
      vec support = regspace<vec>(1, n_items);

      // Vector which holds the proposed ranks
      vec ranks = zeros(n_items);
      vec ranks2 = zeros(n_items);

      // Probability of the ranks we get
      double log_q = 0;

      // n_items random uniform numbers
      vec u(n_items);
      for(int i{}; i < n_items; i++) u(i) = std::log(R::runif(0, 1));

      // Loop over possible values given to item j in random order
      Rcpp::IntegerVector a = Rcpp::seq(0, n_items - 1);
      a = Rcpp::sample(a, n_items);
      ivec myind = Rcpp::as<ivec>(Rcpp::wrap(a));

      for(int j = 0; j < n_items; ++j){
        int jj = myind(j);
        // Find the elements that have not been taken yet
        uvec inds = find(support != 0);
        vec log_prob(inds.n_elem);

        // Number of elements
        int k_max = inds.n_elem;

        // Reference vector
        vec r1 = inds + ones(k_max);
        // Sampled vector
        vec r2 = rho(jj) * ones(k_max);
        // Probability of sample. Note that this is a vector quantity.
        log_prob = - alpha / n_items * pow(abs(r1 - r2), (metric == "footrule") ? 1. : 2.);
        log_prob = log_prob - std::log(accu(exp(log_prob)));
        vec log_cpd = log(cumsum(exp(log_prob)));

        // Draw a random sample
        int item_index = as_scalar(find(log_cpd > u(jj), 1));
        ranks(jj) = as_scalar(inds(item_index)) + 1;

        log_q += log_prob(item_index);

        // Remove the realized rank from support
        support(ranks(jj) - 1) = 0;
      }

      // Increment the partition function
      partfun(i) = - alpha / n_items * distfun->d(ranks, rho) - log_q;
    }
    // Average over the Monte Carlo samples
    // Using this trick: https://www.xarg.org/2016/06/the-log-sum-exp-trick-in-machine-learning/
    double maxval = max(partfun);
    logZ(t) = maxval + std::log(accu(exp(partfun - maxval))) - std::log(static_cast<double>(nmc));
  }
  return logZ;
}
