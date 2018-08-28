#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::vec find_pairwise_limits(int left_limit, int right_limit, int element,
                               arma::uvec ordering, arma::vec possible_rankings){


  // Find the indices of the constrained elements which are preferred to *element*
  // Find the index of this element
  arma::uvec element_ind = arma::find(ordering == element, 1);

  // If there are any elements above, we must find the possible rankings
  if(element_ind(0) < (ordering.n_elem - 1)){

    arma::uvec preferred_element_inds = arma::regspace<arma::uvec>(element_ind(0) + 1, 1, ordering.n_elem - 1);
    arma::uvec preferred_elements = ordering(preferred_element_inds);

    arma::vec rankings_above = possible_rankings.elem(preferred_elements - 1);

    left_limit = arma::max(rankings_above);

  }

  // If there are any elements below, we must find the possible rankings
  if(element_ind(0) > 0){
    arma::uvec disfavored_element_inds = arma::regspace<arma::uvec>(0, 1, element_ind(0) - 1);

    arma::uvec disfavored_elements = ordering(disfavored_element_inds);
    arma::vec rankings_below = possible_rankings.elem(disfavored_elements - 1);

    right_limit = arma::min(rankings_below);
  }

  arma::vec result(2);
  result(0) = left_limit;
  result(1) = right_limit;
  return result;
}

// [[Rcpp::export]]
arma::vec find_pairwise_limits_old(int left_limit, int right_limit, int element,
                                   arma::mat pairwise_preferences,
                                   arma::vec possible_rankings, int i){
  // Find the indices of the constrained elements which are preferred to *element*
  arma::uvec preferred_element_inds = arma::intersect(
    arma::find(pairwise_preferences.col(0) == (i + 1)),
    arma::find(pairwise_preferences.col(1) == element) // bottom_item is element
  );

  arma::uvec preferred_elements = arma::conv_to<arma::uvec>::from(pairwise_preferences.col(2));
  preferred_elements = preferred_elements.elem(preferred_element_inds);

  arma::vec rankings_above = possible_rankings.elem(preferred_elements - 1);

  if(possible_rankings.n_elem > 0) {
    left_limit = arma::max(possible_rankings);
  }

  // Find the indices of the constrained elements which are disfavored to *element*
  arma::uvec disfavored_element_inds = arma::intersect(
    arma::find(pairwise_preferences.col(0) == (i + 1)),
    arma::find(pairwise_preferences.col(2) == element) // top_item is element
  );

  arma::uvec disfavored_elements = arma::conv_to<arma::uvec>::from(pairwise_preferences.col(1));
  disfavored_elements = disfavored_elements.elem(disfavored_element_inds);

  arma::vec rankings_below = possible_rankings.elem(disfavored_elements - 1);

  if(possible_rankings.n_elem > 0) {
    right_limit = arma::min(possible_rankings);
  }

  arma::vec result(2);
  result(0) = left_limit;
  result(1) = right_limit;
  return result;
}

void augment_pairwise(
    arma::mat& rankings,
    const arma::umat& cluster_assignment,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& linear_ordering,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::vec& aug_acceptance,
    const bool& clustering,
    bool& augmentation_accepted
){

  for(int i = 0; i < n_assessors; ++i){

    // Draw an integer between 1 and n_items
    int element = arma::randi<int>(arma::distr_param(1, n_items));

    arma::uvec ordering = linear_ordering[i];

    bool element_is_constrained = arma::any(ordering == element);

    // Left and right limits of the interval we draw ranks from
    // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
    int left_limit = 0, right_limit = n_items + 1;

    if(element_is_constrained){

      arma::vec limits = find_pairwise_limits(left_limit, right_limit, element, ordering, rankings.col(i));
      left_limit = limits(0);
      right_limit = limits(1);
    }

    // Now complete the leap step by drawing a new proposal uniformly between
    // right_limit + 1 and left_limit - 1
    int proposed_rank = arma::randi<int>(arma::distr_param(left_limit + 1, right_limit - 1));

    // Assign the proposal to the (element-1)th element
    arma::vec proposal = rankings.col(i);
    proposal(element - 1) = proposed_rank;

    double delta_r;
    arma::uvec indices;

    // Do the shift step
    shift_step(proposal, rankings.col(i), element, delta_r, indices);

    // Finally, decide whether to accept the proposal or not
    // Draw a uniform random number
    double u = log(arma::randu<double>());

    // Find which cluster the assessor belongs to
    int cluster;
    if(clustering){
      cluster = cluster_assignment(i, t);
    } else {
      cluster = 0;
    }

    double ratio = -alpha(cluster) / n_items *
      (get_rank_distance(proposal, rho.col(cluster), metric) -
      get_rank_distance(rankings.col(i), rho.col(cluster), metric));

    if(ratio > u){
      rankings.col(i) = proposal;
      ++aug_acceptance(i);
      augmentation_accepted = true;
    } else {
      augmentation_accepted = false;
    }
  }

}

