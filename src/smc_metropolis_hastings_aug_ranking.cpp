// #include <RcppArmadillo.h>
// #include "smc.h"
// #include "misc.h"
// #include "setdiff.h"
//
// using namespace arma;
//
// // [[Rcpp::depends(RcppArmadillo)]]
//
// arma::vec metropolis_hastings_aug_ranking(
// 	const double& alpha,
// 	const arma::vec& rho,
// 	const int& n_items,
// 	const arma::vec& partial_ranking,
// 	const arma::vec& current_ranking,
// 	const bool& pseudo,
// 	const uvec& missing_indicator,
// 	const std::string& metric = "footrule"
// ) {
//   double forward_backward_prob{};
//   vec proposed_ranking;
//
//   // find items missing from original observed ranking
//
//   uvec missing_ranks = find(missing_indicator == 1);
//
//   // find unallocated ranks from original observed ranking
//   const vec remaining_set = partial_ranking(missing_ranks);
//
//   if(pseudo){
//     // randomly permute the unranked items to give the order in which they will
//     // be allocated
//     uvec item_ordering = shuffle(missing_ranks);
//
//     // Calculate probabilities
//     const Rcpp::List proposal = calculate_forward_probability(
//       current_ranking,
//       item_ordering,
//       missing_indicator,
//       alpha, rho, metric
//     );
//     const vec& proposed_augmented_ranking = proposal["aug_ranking"];
//     double forward_prob = proposal["forward_prob"];
//
//     double backward_prob = calculate_backward_probability(
//       item_ordering, partial_ranking, current_ranking, remaining_set, rho,
//       alpha, n_items, metric
//     );
//     forward_backward_prob = - std::log(forward_prob) + std::log(backward_prob);
//     proposed_ranking = proposed_augmented_ranking;
//
//   } else {
//     // Subset by element position and set equal to the now permuted remaining set
//     proposed_ranking = partial_ranking;
//     proposed_ranking.elem(missing_ranks) = shuffle(remaining_set);
//
//     // set the augmented partial ranking as the proposed augmented ranking
//   }
//
//
//   /* MH TIME ------------------------------------------------------------- */
//   // Calculate the log posterior of the current and proposed rankings.
//   // NB the current can usually be stored to save recalculating it, but we're not caring about that yet
//   const double curr_logpost = get_exponent_sum(     \
//     alpha, rho, n_items, current_ranking, metric\
//   );
//   const double prop_logpost = get_exponent_sum(                \
//     alpha, rho, n_items, proposed_ranking, metric          \
//   );
//
//   const double log_acceptance_prob = prop_logpost - curr_logpost + forward_backward_prob;
//
//   if (std::log(Rcpp::as<double>(Rcpp::runif(1, 0, 1))) < log_acceptance_prob) {
//     return(proposed_ranking);
//   } else {
//     return(current_ranking);
//   }
// }
