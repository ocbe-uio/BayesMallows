#include "RcppArmadillo.h"
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Get Sample Probabilities
//' @description Calculate probability of assigning a set of specific ranks to an specific item
//' given its rank in the consensus ranking
//'
//' @param rho_item_rank An integer value rank of an item in the current consensus ranking
//' @param alpha Numeric value og the scale parameter
//' @param remaining_set_rank A sequence of integer values of the set of possible ranks that we can assign the item
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @param n_items Integer is the number of items in the consensus ranking
//' @return sample_prob_list A numeric sequence of sample probabilities for selecting a specific rank given the current
//'         rho_item_rank
//' @export
// get_sample_probabilities = function(rho_item_rank, alpha, remaining_set_ranks, metric, n_items){


//   # define a set of probs list
//   num_ranks = length(remaining_set_ranks)
//   sample_prob_list = rep(0, num_ranks)

//   # cycle through each item and calculate its specific prob
//   for (ii in 1:num_ranks){
//     item_rank = remaining_set_ranks[ii]
//     sample_prob = (-(alpha/n_items)*BayesMallows:::get_rank_distance(rho_item_rank, item_rank, metric) )
//     sample_prob_list[ii] = sample_prob
//   }

//   # normalise probs
//   maxw = max(sample_prob_list)
//   w = exp(sample_prob_list-maxw)
//   sample_prob_list = w/sum(w)

//   return(sample_prob_list)
// }
