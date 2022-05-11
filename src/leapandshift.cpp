#include <RcppArmadillo.h>
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
                    const vec& rho, int leap_size, bool reduce_indices){

  // Set proposal equal to current
  rho_proposal = rho;

  // Help vectors
  vec support;

  // Number of items
  int n = rho.n_elem;

  // Other helper variables
  double support_new;

  // Leap 1
  // 1, sample u randomly between 1 and n
  int u = randi<int>(distr_param(0, n - 1));


  // 2, compute the set S for sampling the new rank
  double dobL = static_cast<double>(leap_size);
  double dobn = static_cast<double>(n);

  // Defining linspace lengths here to avoid duplication in code
  double length1 = std::min(rho(u) - 1, dobL);
  double length2 = std::min(n - rho(u), dobL);

  if ((rho(u) > 1) && (rho(u) < n)) {
    support = join_cols(
      linspace(std::max(1.0, rho(u) - leap_size), rho(u) - 1, length1),
      linspace(rho(u) + 1, std::min(dobn, rho(u) + leap_size), length2));
  } else if(rho(u) == 1){
    support = linspace(rho(u) + 1, std::min(dobn, rho(u) + leap_size), length2);
  } else if(rho(u) == n){
    support = linspace(std::max(1.0, rho(u) - leap_size), rho(u) - 1, length1);
  }

  // 3. assign a random element of the support set, this completes the leap step
  int index = randi<int>(distr_param(0, support.n_elem-1));
  // Picked element index-1 from the support set
  rho_proposal(u) = support(index);



  // Compute the associated transition probabilities
  if(std::abs(rho_proposal(u) - rho(u)) == 1){
    // in this case the transition probabilities coincide! (and in fact for leap_size = 1 the L&S is symmetric)
    support_new = std::min(rho_proposal(u) - 1, dobL) + std::min(n - rho_proposal(u), dobL);
    prob_forward = 1.0 / (n * support.n_elem) + 1.0 / (n * support_new);
    prob_backward = prob_forward;
  } else {
    // P(proposed|current)
    prob_forward = 1.0 / (n * support.n_elem);
    // P(current|proposed)
    support_new = std::min(rho_proposal(u) - 1, dobL) + std::min(n - rho_proposal(u), dobL);
    prob_backward = 1.0 / (n * support_new);
  }


  shift_step(rho_proposal, rho, u, indices);

  if(!reduce_indices){
    indices = regspace<uvec>(0, n - 1);
  }
}
