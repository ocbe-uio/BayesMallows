#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]

void augment_pairwise(
    arma::mat& R,
    const double& alpha,
    const arma::vec& rho,
    const std::string& metric,
    const arma::mat& pairwise_preferences,
    const arma::mat& constrained_elements,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::mat& aug_acceptance,
    int& aug_diag_index,
    const int& aug_diag_thinning
){

  for(int i = 0; i < n_assessors; ++i){

    // Draw an integer between 1 and n_items
    int element = arma::as_scalar(arma::randi(1, arma::distr_param(1, n_items)));

    // Find the constrained elements for this particular assessor
    // First, check which columns the assessor has
    // i + 1 because the assessor numbering starts at 1
    arma::uvec inds = arma::find(constrained_elements.row(0) == (i + 1));

    // Then find the corresponding elements
    arma::vec constr = constrained_elements.row(1).t();
    constr = constr.elem(inds);

    // Check if the drawn element is in the constraint set
    arma::uvec element_ind = arma::find(constr == element);

    // Lower and upper limits of the rankings we will draw
    int lower_limit = 0, upper_limit = n_items + 1;

    if(element_ind.n_elem > 0){

      // We extract the pairwise preferences of assessor i
      arma::uvec pref_inds = arma::find(pairwise_preferences.row(0) == (i + 1));
      arma::mat prefs = pairwise_preferences.submat(1, arma::min(pref_inds),
                                                    2, arma::max(pref_inds));

      // We must find the items that are ranked below the drawn element
      // First, where is this element preferred
      arma::uvec rank_inds = arma::find(prefs.row(1) == element);
      // Next, which elements are ranked below
      arma::vec elements_below = prefs.row(0).t();
      elements_below = elements_below.elem(rank_inds);
      // Subtract 1 to go from items to indices
      arma::uvec indices_below = arma::conv_to<arma::uvec>::from(elements_below) -
        arma::ones<arma::uvec>(elements_below.n_elem);

      // Extract the ranks of the elements below
      arma::vec ranks_below;
      if(rank_inds.n_elem > 0){
        ranks_below = R.col(i);
        ranks_below = ranks_below.elem(indices_below);
        lower_limit = arma::max(ranks_below);
      }

      // Next, we find the items that are ranked above the drawn element
      // First, where is this element not preferred
      rank_inds = arma::find(prefs.row(0) == element);
      // Next, which elements are ranked above
      arma::vec elements_above = prefs.row(1).t();
      elements_above = elements_above.elem(rank_inds);
      // Subtract 1 to go from items to indices
      arma::uvec indices_above = arma::conv_to<arma::uvec>::from(elements_above) -
        arma::ones<arma::uvec>(elements_above.n_elem);

      // Extract the ranks of the elements above
      arma::vec ranks_above;
      if(rank_inds.n_elem > 0){
        ranks_above = R.col(i);
        ranks_above = ranks_above.elem(indices_above);
        upper_limit = arma::min(ranks_above);
      }

    }

    // Now complete the leap step by drawing a new proposal uniformly between
    // lower_limit + 1 and upper_limit - 1

    int proposed_element = arma::as_scalar(
      arma::randi(1, arma::distr_param(lower_limit + 1, upper_limit - 1)));


    // Assign the proposal to the (element-1)th element
    arma::vec proposal = R.col(i);
    proposal(element - 1) = proposed_element;

    double delta_r;
    arma::vec indices;

    // Do the shift step
    shift_step(proposal, R.col(i), element, delta_r, indices);

    // Finally, decide whether to accept the proposal or not
    // Draw a uniform random number
    double u = log(arma::randu<double>());

    double ratio = -alpha / n_items *
      (get_rank_distance(proposal, rho, metric) -
      get_rank_distance(R.col(i), rho, metric));

    if(ratio > u){
      R.col(i) = proposal;
      ++aug_acceptance(i, aug_diag_index);
    }
  }

  // If appropriate, increment the augmentation diagnostic index
  if(t % aug_diag_thinning == 0){
    ++aug_diag_index;
  }
}
