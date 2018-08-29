#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]

void find_pairwise_limits(int& left_limit, int& right_limit, const int& element,
                          const arma::uvec& ordering,
                          const arma::vec& current_ranking){

  // Find the indices of the constrained elements which are preferred to *element*
  // Find the index of this element
  arma::uvec element_ind = arma::find(ordering == element, 1);

  // If there are any elements above, we must find the possible rankings
  if(element_ind(0) < (ordering.n_elem - 1)){

    arma::uvec preferred_element_inds = arma::regspace<arma::uvec>(element_ind(0) + 1, 1, ordering.n_elem - 1);
    arma::uvec preferred_elements = ordering(preferred_element_inds);

    arma::vec rankings_above = current_ranking.elem(preferred_elements - 1);

    left_limit = arma::max(rankings_above);

  }

  // If there are any elements below, we must find the possible rankings
  if(element_ind(0) > 0){
    arma::uvec disfavored_element_inds = arma::regspace<arma::uvec>(0, 1, element_ind(0) - 1);

    arma::uvec disfavored_elements = ordering(disfavored_element_inds);
    arma::vec rankings_below = current_ranking.elem(disfavored_elements - 1);

    right_limit = arma::min(rankings_below);
  }

}

//' Find the left and right limit for proposing augmented rank.
//'
//' This function implements a part of the modified leap-and-shift
//' algorithm for proposing augmented ranks that agree with the
//' transitive closure of pairwise preferences stated by the assessor. The
//' algorithm is described on pp. 21-22 of \insertCite{vitelli2018;textual}{BayesMallows}.
//' \code{find_pairwise_limits} returns the left and right limits, \eqn{l_{j}} and
//' \eqn{r_{j}}. The proposed new ranking of the given element is then
//' sampled uniformly from the set \eqn{\{l_{j} + 1, \dots, r_{j} - 1\}}, whereupon
//' a shift step is applied. This function is separated out mainly for testing
//' purposes.
//'
//' @param u An integer specifying the element to modify. In the modified leap-and-shift
//' algorithm, \eqn{u} is drawn randomly from the set \eqn{\{1, \dots, n_{items}\}}.
//'
//' @param ordering A vector that contains the linear ordering of elements implied by
//' the preferences stated by the assessor. Note that the unconstrained elements are
//' not supposed to be part of this vector.
//'
//' @param current_ranking A vector that contains the current complete ranking for
//' the assessor. This correspondings to \eqn{\tilde{R}_{j}} in the modified leap-and-shift
//' algorithm.
//'
//'
//' @return
//' A vector of size 2. The first element is the left limit (\eqn{l_{j}}) and
//' the second element is the right limit (\eqn{r_{j}}).
//'
//' @references \insertAllCited{}
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec find_pairwise_limits(int u, arma::uvec ordering, arma::vec current_ranking){
  int n_items = current_ranking.n_elem;

  int left_limit = 0;
  int right_limit = n_items + 1;

  bool element_is_constrained = arma::any(ordering == u);
  if(element_is_constrained){
    find_pairwise_limits(left_limit, right_limit, u, ordering, current_ranking);
  }


  arma::vec result(2); result(0) = left_limit; result(1) = right_limit;
  return result;
}



void propose_pairwise_augmentation(arma::vec& proposal,
                                   const arma::mat& rankings,
                                   const Rcpp::List& linear_ordering,
                                   const int& n_items,
                                   const int& i){
  // Draw an integer between 1 and n_items
  int element = arma::randi<int>(arma::distr_param(1, n_items));

  arma::uvec ordering = linear_ordering[i];

  bool element_is_constrained = arma::any(ordering == element);

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;

  if(element_is_constrained){
    find_pairwise_limits(left_limit, right_limit, element,
                         ordering, rankings.col(i));
  }

  // Now complete the leap step by drawing a new proposal uniformly between
  // right_limit + 1 and left_limit - 1
  int proposed_rank = arma::randi<int>(arma::distr_param(left_limit + 1, right_limit - 1));

  // Assign the proposal to the (element-1)th element
  proposal = rankings.col(i);
  proposal(element - 1) = proposed_rank;

  double delta_r;
  arma::uvec indices;

  // Do the shift step
  shift_step(proposal, rankings.col(i), element, delta_r, indices);

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

    arma::vec proposal;
    propose_pairwise_augmentation(proposal, rankings, linear_ordering, n_items, i);

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

// [[Rcpp::export]]
arma::mat check_pairwise_augmentation(arma::mat& rankings,
                                      Rcpp::List linear_ordering){

  arma::vec proposal;
  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    propose_pairwise_augmentation(proposal, rankings, linear_ordering,
                                  n_items, i);
    rankings.col(i) = proposal;
  }

  return rankings;
}

