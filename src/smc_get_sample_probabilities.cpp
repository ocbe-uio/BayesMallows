#include "RcppArmadillo.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Get Sample Probabilities
//' @description Calculate probability of assigning a set of specific ranks to an specific item
//' given its rank in the consensus ranking
//'
//' @param rho_item_rank An integer value rank of an item in the current consensus ranking
//' @param alpha Numeric value of the scale parameter
//' @param remaining_set_ranks A sequence of integer values of the set of possible ranks that we can assign the item
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @param n_items Integer is the number of items in the consensus ranking
//' @return sample_prob_list A numeric sequence of sample probabilities for selecting a specific rank given the current
//'         rho_item_rank
//' @export
// [[Rcpp::export]]
arma::vec get_sample_probabilities(
  const arma::vec rho_item_rank,
  const double alpha,
  const arma::vec remaining_set_ranks,
  const std::string metric,
  const int n_items
) {
  // define a set of probs list
  const unsigned int& num_ranks = remaining_set_ranks.n_elem;
  vec sample_prob_list = Rcpp::rep(0.0, num_ranks);

  // cycle through each item and calculate its specific prob
  for (uword ii = 0; ii < num_ranks; ++ii) {
    const vec item_rank = {remaining_set_ranks(ii)};
    const double rank_dist = get_rank_distance(rho_item_rank, item_rank, metric);
    const double sample_prob = -(alpha / n_items) * rank_dist;
    sample_prob_list(ii) = sample_prob;
  }

  // normalise probs
  const double maxw = max(sample_prob_list);
  const vec w = exp(sample_prob_list - maxw);
  sample_prob_list = w / arma::sum(w);

  return(sample_prob_list);
}
