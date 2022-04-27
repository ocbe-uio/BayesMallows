#include "RcppArmadillo.h"
#include "smc.h"
#include "misc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Calculate Backward Probability
//' @description Function to calculate probability of assigning a set of specific ranks to an specific item
//' given its rank in the consensus ranking
//'
//' @param item_ordering A vector of integer values to represent the specified queue of which unranked item to assign a rank for the proposed augmented ranking
//' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obtained from calculate_forward_probability function
//' @param remaining_set A vector of integer values to represent the elements (ranks) missing from original observed ranking
//' @param rho Numeric vector specifying the consensus ranking
//' @param alpha Numeric value of the scale parameter
//' @param n_items Integer is the number of items in a ranking
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return backward_auxiliary_ranking_probability A numerical value of creating the previous augmented ranking using the same item ordering used to create the
//' new augmented ranking in calculate_forward_probability function.
//' @export
// [[Rcpp::export]]
double calculate_backward_probability(
  arma::uvec item_ordering,
  arma::vec partial_ranking,
  arma::vec current_ranking,
  arma::vec remaining_set,
  const arma::vec rho,
  const double alpha,
  const int n_items,
  const std::string metric
) {
  // given an old and new item ordering, sample ranking with new ordering and
  // calc the forward and backward prob

  // show the augmented parts of the current ranking
  // item ordering is the order of which items are assigned ranks in a specified
  // order
  const uword& num_items_unranked = item_ordering.n_elem;

  // initialise prob of creating augmented ranking
  double backward_auxiliary_ranking_probability = 1.0;

  // Adjust item_ordering depending on whether it uses R or C++ indices
  item_ordering = maybe_offset_indices(partial_ranking, item_ordering);

  if (num_items_unranked != 1) {
    // show the augmented parts of the current ranking
    current_ranking = current_ranking.elem(item_ordering);

  /* ====================================================== */
  /* LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY     */
  /* ====================================================== */
  // given the old and new item ordering and the list of missing rank, determine
  // the sample probs for each iteration
    for (uword jj = 0; jj < num_items_unranked - 1; ++jj) {

      // items to sample rank
      const uword& item_to_sample_rank = item_ordering(jj);

      // the rank of item in rho
      vec rho_item_rank;
      rho_item_rank = rho(item_to_sample_rank);

      // next we get the sample probabilites for selecting a particular rank for
      // an item based on the current alpha and the rho rank for that item
      const vec& sample_prob_list = get_sample_probabilities(\
        rho_item_rank, alpha, remaining_set, metric, n_items\
      );

      // save the probability of selecting the specific item rank in the old
      // augmented ranking
      const uvec& sample_prob = find(remaining_set == current_ranking(jj));
      backward_auxiliary_ranking_probability = \
        backward_auxiliary_ranking_probability * \
        arma::as_scalar(sample_prob_list(sample_prob));

      // remove selected auxiliary rank from the set of remaining possibles
      // ranks to select
      remaining_set = remaining_set(find(remaining_set != current_ranking(jj)));
    }
  }
  return backward_auxiliary_ranking_probability;
}
