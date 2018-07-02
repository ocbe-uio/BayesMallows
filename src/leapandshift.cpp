#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


Rcpp::List leap_and_shift(arma::vec rho, int L){

  // Declare the proposed rank vector
  arma::vec proposal = rho;

  // Help vectors
  arma::vec support, indices;

  // Number of items
  int n = rho.n_elem;

  // Other helper variables
  int u, index;
  double delta_r, prob_forward, prob_backward, support_new;

  // Leap 1
  // 1, sample u randomly between 1 and n
  u = arma::as_scalar(arma::randi(1, arma::distr_param(1, n)));

  // 2, compute the set S for sampling the new rank
  // Defining versions of L and n converted to double, to avoid duplication in code
  double dobL = static_cast<double>(L);
  double dobn = static_cast<double>(n);

  // Defining linspace lengths here to avoid duplication in code
  double length1 = std::min(rho(u - 1) - 1, dobL);
  double length2 = std::min(n - rho(u - 1), dobL);

  if((rho(u - 1) > 1) & (rho(u - 1) < n)){
    support = arma::join_cols(
      arma::linspace(
        std::max(1.0, rho(u - 1) - L), rho(u - 1) - 1, length1
      ),
      arma::linspace(
        rho(u - 1) + 1, std::min(dobn, rho(u - 1) + L), length2
      )
    );
  } else if(rho(u - 1) == 1){
    support = arma::linspace(
      rho(u - 1) + 1,
      std::min(dobn, rho(u - 1) + L),
      length2
    );

  } else if(rho(u - 1) == n){
    support = arma::linspace(
      std::max(1.0, rho(u - 1) - L), rho(u - 1) - 1,
      length1);
  }

  // 3. assign a random element of the support set, this completes the leap step
  index = arma::as_scalar(arma::randi(1, arma::distr_param(0, support.n_elem-1)));
  // Picked element index-1 from the support set
  proposal(u-1) = support(index);

  // Compute the associated transition probabilities (BEFORE THE SHIFT STEP, WHICH IS DETERMINISTIC --> EASIER)
  if(std::abs(proposal(u - 1) - rho(u - 1)) == 1){
    // in this case the transition probabilities coincide! (and in fact for L = 1 the L&S is symmetric)
    support_new = std::min(proposal(u - 1) - 1, dobL) + std::min(n - proposal(u - 1), dobL);
    prob_forward = 1.0 / (n * support.n_elem) + 1.0 / (n * support_new);
    prob_backward = prob_forward;
  } else {
    // P(proposed|current)
    prob_forward = 1.0 / (n * support.n_elem);
    // P(current|proposed)
    support_new = std::min(proposal(u - 1) - 1, dobL) + std::min(n - proposal(u-1), dobL);
    prob_backward = 1.0 / (n * support_new);
  }

  // Shift step:
  delta_r = proposal(u - 1) - rho(u - 1);
  indices = arma::zeros(std::abs(delta_r) + 1);
  indices[0] = u-1;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k =- 1; k>=delta_r; --k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      proposal(index) += 1;
      indices[-(k)] = index;
    }
  }


  return Rcpp::List::create(Rcpp::Named("proposal") = proposal,
                            Rcpp::Named("indices") = indices,
                            Rcpp::Named("delta_r") = delta_r,
                            Rcpp::Named("prob_forward") = prob_forward,
                            Rcpp::Named("prob_backward") = prob_backward
  );
}


