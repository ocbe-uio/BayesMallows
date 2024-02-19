#include <RcppArmadillo.h>
#include "distances.h"
#include "leapandshift.h"
using namespace arma;

void LeapShiftObject::shift_step(const vec& rho, int u){
  // Shift step:
  double delta_r = rho_proposal(u) - rho(u);
  indices = zeros<uvec>(std::abs(delta_r) + 1);
  indices[0] = u;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      int index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      int index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) += 1;
      indices[-(k)] = index;
    }
  }
}

LeapShiftObject leap_and_shift(
    const vec& rho, int leap_size, const std::unique_ptr<Distance>& distfun){
  LeapShiftObject ls{rho};
  int n_items = rho.n_elem;

  // Leap 1
  // 1, sample u randomly between 1 and n_items
  ivec a = Rcpp::sample(n_items, 1) - 1;
  int u = a(0);

  // 2, compute the set S for sampling the new rank
  vec support = join_cols(regspace(std::max(1.0, rho(u) - leap_size), 1, rho(u) - 1),
    regspace(rho(u) + 1, 1, std::min(n_items * 1.0, rho(u) + leap_size)));

  // 3. assign a random element of the support set, this completes the leap step
  ivec b = Rcpp::sample(support.n_elem, 1) - 1;
  int index = b(0);
  // Picked element index-1 from the support set
  ls.rho_proposal(u) = support(index);

  double support_new = std::min(ls.rho_proposal(u) - 1, leap_size * 1.0) +
    std::min(n_items - ls.rho_proposal(u), leap_size * 1.0);

  ls.shift_step(rho, u);
  distfun->update_leap_and_shift_indices(ls.indices, n_items);

  if(std::abs(ls.rho_proposal(u) - rho(u)) == 1){
    ls.prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new);
    ls.prob_backward = ls.prob_forward;
  } else {
    ls.prob_forward = 1.0 / (n_items * support.n_elem);
    ls.prob_backward = 1.0 / (n_items * support_new);
  }
  return ls;
}
