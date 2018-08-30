#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]

void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const arma::vec& current_ranking){

  // Find the items which are preferred to the given item
  // Items go from 1, ..., n_items, so must use [item - 1]
  arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[item - 1]);
  arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[item - 1]);

  // If there are any items above, we must find the possible rankings
  if(items_above.n_elem > 0){
    // Again subtracting 1 because of zero-first indexing
    // Find all the rankings of the items that are preferred to *item*
    arma::vec rankings_above = current_ranking.elem(items_above - 1);
    left_limit = arma::max(rankings_above);
  }

  // If there are any items below, we must find the possible rankings
  if(items_below.n_elem > 0){
    // Find all the rankings of the items that are disfavored to *item*
    arma::vec rankings_below = current_ranking.elem(items_below - 1);
    right_limit = arma::min(rankings_below);
  }

}



void propose_pairwise_augmentation(arma::vec& proposal,
                                   const arma::mat& rankings,
                                   const Rcpp::List& constraints,
                                   const int& n_items,
                                   const int& i){

  // Extract the constraints for this particular assessor
  Rcpp::List assessor_constraints = Rcpp::as<Rcpp::List>(constraints[i]);
  arma::uvec constrained_items = Rcpp::as<arma::uvec>(assessor_constraints[0]);

  // Draw an integer between 1 and n_items
  int item = arma::randi<int>(arma::distr_param(1, n_items));
  // Check if the item is constrained for this assessor
  bool item_is_constrained = arma::any(constrained_items == item);

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;

  if(item_is_constrained){
    find_pairwise_limits(left_limit, right_limit, item,
                         assessor_constraints, rankings.col(i));
  }

  // Now complete the leap step by drawing a new proposal uniformly between
  // left_limit + 1 and right_limit - 1
  int proposed_rank = arma::randi<int>(arma::distr_param(left_limit + 1, right_limit - 1));

  // Assign the proposal to the (item-1)th item
  proposal = rankings.col(i);
  proposal(item - 1) = proposed_rank;

  double delta_r;
  arma::uvec indices;

  // Do the shift step
  shift_step(proposal, rankings.col(i), item, delta_r, indices);

}

void augment_pairwise(
    arma::mat& rankings,
    const arma::umat& cluster_assignment,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::vec& aug_acceptance,
    const bool& clustering,
    bool& augmentation_accepted
){

  for(int i = 0; i < n_assessors; ++i){
    // Call the function which creates a proposal
    arma::vec proposal;
    propose_pairwise_augmentation(proposal, rankings, constraints, n_items, i);

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

