#include <RcppArmadillo.h>
#include "distances.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void shift_step(vec& rho_proposal, const vec& rho,
                const int& u, uvec& indices){
  // Shift step:
  double delta_r = rho_proposal(u) - rho(u);
  indices = zeros<uvec>(std::abs(delta_r) + 1);
  indices[0] = u;
  int index;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) += 1;
      indices[-(k)] = index;
    }
  }
}


void leap_and_shift(vec& rho_proposal, uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const vec& rho, int leap_size,
                    const std::unique_ptr<Distance>& distfun){
  rho_proposal = rho;
  vec support;
  int n_items = rho.n_elem;
  double support_new;

  // Leap 1
  // 1, sample u randomly between 1 and n_items
  ivec a = Rcpp::sample(n_items, 1) - 1;
  int u = a(0);

  // 2, compute the set S for sampling the new rank
  support = join_cols(regspace(std::max(1.0, rho(u) - leap_size), 1, rho(u) - 1),
    regspace(rho(u) + 1, 1, std::min(n_items * 1.0, rho(u) + leap_size)));

  // 3. assign a random element of the support set, this completes the leap step
  ivec b = Rcpp::sample(support.n_elem, 1) - 1;
  int index = b(0);
  // Picked element index-1 from the support set
  rho_proposal(u) = support(index);

  support_new = std::min(rho_proposal(u) - 1, leap_size * 1.0) + std::min(n_items - rho_proposal(u), leap_size * 1.0);
  if(std::abs(rho_proposal(u) - rho(u)) == 1){
    prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new);
    prob_backward = prob_forward;
  } else {
    prob_forward = 1.0 / (n_items * support.n_elem);
    prob_backward = 1.0 / (n_items * support_new);
  }
  shift_step(rho_proposal, rho, u, indices);
  distfun->update_leap_and_shift_indices(indices, n_items);
}
