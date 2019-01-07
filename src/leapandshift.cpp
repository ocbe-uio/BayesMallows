#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

void shift_step(arma::vec& rho_proposal, const arma::vec& rho,
                const int& u, double& delta_r, arma::uvec& indices){
  // Shift step:
  delta_r = rho_proposal(u - 1) - rho(u - 1);
  indices = arma::zeros<arma::uvec>(std::abs(delta_r) + 1);
  indices[0] = u-1;
  int index;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      rho_proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      rho_proposal(index) += 1;
      indices[-(k)] = index;
    }
  }
}


void leap_and_shift(arma::vec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::vec& rho, int leap_size, bool reduce_indices){

  // Set proposal equal to current
  rho_proposal = rho;

  // Help vectors
  arma::vec support;

  // Number of items
  int n = rho.n_elem;

  // Other helper variables
  int u, index;
  double delta_r, support_new;

  // Leap 1
  // 1, sample u randomly between 1 and n
  u = arma::as_scalar(arma::randi(1, arma::distr_param(1, n)));

  // 2, compute the set S for sampling the new rank
  double dobL = static_cast<double>(leap_size);
  double dobn = static_cast<double>(n);

  // Defining linspace lengths here to avoid duplication in code
  double length1 = std::min(rho(u - 1) - 1, dobL);
  double length2 = std::min(n - rho(u - 1), dobL);

  if((rho(u - 1) > 1) & (rho(u - 1) < n)){
    support = arma::join_cols(arma::linspace(
      std::max(1.0, rho(u - 1) - leap_size), rho(u - 1) - 1, length1),
      arma::linspace(rho(u - 1) + 1, std::min(dobn, rho(u - 1) + leap_size), length2));
  } else if(rho(u - 1) == 1){
    support = arma::linspace(rho(u - 1) + 1, std::min(dobn, rho(u - 1) + leap_size), length2);
  } else if(rho(u - 1) == n){
    support = arma::linspace(std::max(1.0, rho(u - 1) - leap_size), rho(u - 1) - 1, length1);
  }

  // 3. assign a random element of the support set, this completes the leap step
  index = arma::as_scalar(arma::randi(1, arma::distr_param(0, support.n_elem-1)));
  // Picked element index-1 from the support set
  rho_proposal(u-1) = support(index);

  // Compute the associated transition probabilities
  if(std::abs(rho_proposal(u - 1) - rho(u - 1)) == 1){
    // in this case the transition probabilities coincide! (and in fact for leap_size = 1 the L&S is symmetric)
    support_new = std::min(rho_proposal(u - 1) - 1, dobL) + std::min(n - rho_proposal(u - 1), dobL);
    prob_forward = 1.0 / (n * support.n_elem) + 1.0 / (n * support_new);
    prob_backward = prob_forward;
  } else {
    // P(proposed|current)
    prob_forward = 1.0 / (n * support.n_elem);
    // P(current|proposed)
    support_new = std::min(rho_proposal(u - 1) - 1, dobL) + std::min(n - rho_proposal(u-1), dobL);
    prob_backward = 1.0 / (n * support_new);
  }

  shift_step(rho_proposal, rho, u, delta_r, indices);

  if(!reduce_indices){
    indices = arma::regspace<arma::uvec>(0, n - 1);
  }
}



