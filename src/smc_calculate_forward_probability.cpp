#include "RcppArmadillo.h"
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Calculate Forward Probability
//' @description Function to calculate probability of assigning a set of specific ranks to an specific item
//' given its rank in the consensus ranking
//' @export
//'
//' @param item_ordering A vector of integer values to represent the specified queue of which unranked item to assign a rank for the proposed augmented ranking
//' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param remaining_set A vector of integer values to represent the elements (ranks) missing from original observed ranking
//' @param rho Numeric vector specifying the consensus ranking
//' @param alpha Numeric value og the scale parameter
//' @param n_items Integer is the number of items in a ranking
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return List containing aug_ranking, a ranking sequence vector of the proposed augmented ranking and
//'         forward_prob a numerical value of the probability of creating the augmented ranking
//'         using the pseudolikelihood augmentation.
// [[Rcpp::export]]
double calculate_forward_probability(
  arma::vec item_ordering,
  arma::vec partial_ranking,
  arma::vec remaining_set,
  arma::vec rho,
  double alpha,
  int n_items,
  std::string metric
){

  // item ordering is the order of which items are assigned ranks in a specified order
  // num_items_unranked = length(item_ordering)

//   # prob of creating augmented ranking
//   forward_auxiliary_ranking_probability = 1

//   if(num_items_unranked == 1){

//     # create new agumented ranking by sampling remaining ranks from set uniformly
//     partial_ranking[is.na(partial_ranking)] <- remaining_set


//   }else{

//     auxiliary_ranking = rep(0, num_items_unranked)

//     ########################################################
//     ## LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY
//     ########################################################
//     # given the old and new item ordering and the list of missing rank, determine the sample probs for each iteration
//     for (jj in 1:(num_items_unranked-1)){

//       # items to sample rank
//       item_to_sample_rank = item_ordering[jj]

//       # the rank of item in rho
//       rho_item_rank = rho[item_to_sample_rank]

//       # next we get the sample probabilites for selecting a particular rank for an item
//       # based on the current alpha and the rho rank for that item
//       sample_prob_list = get_sample_probabilities(rho_item_rank = rho_item_rank, alpha = alpha,
//                                                   remaining_set_ranks = remaining_set, metric = metric, n_items = n_items)
//       #print(sample_prob_list)

//       # fill in the new augmented ranking going forward
//       auxiliary_ranking[jj] = sample(remaining_set, size = 1, replace = FALSE, prob = sample_prob_list)

//       # save the probability of selecting the specific item rank in the old augmented ranking
//       sample_prob = which(remaining_set == auxiliary_ranking[jj])
//       forward_auxiliary_ranking_probability = forward_auxiliary_ranking_probability * (sample_prob_list[sample_prob])

//       # remove selected auxiliary rank from the set of remaining possibles ranks to select
//       remaining_set = remaining_set[ remaining_set != auxiliary_ranking[jj] ]

//     }

//     # last element in augmented ranking is deterministic - the prob is 1
//     auxiliary_ranking[num_items_unranked] = remaining_set

//     # fit the augmented ranking within the partial rankings with NAs
//     partial_ranking[item_ordering] <- auxiliary_ranking # ranks for items

//   } #end of if else statement


//   output <- list("aug_ranking" = partial_ranking, "forward_prob" =  forward_auxiliary_ranking_probability)
//   return(output)
}