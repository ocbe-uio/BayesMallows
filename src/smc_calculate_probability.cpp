// #include <RcppArmadillo.h>
// #include "smc.h"
// #include "misc.h"
// #include "sample.h"
// #include "missing_data.h"
//
// using namespace arma;
//
// // [[Rcpp::depends(RcppArmadillo)]]
//
//
// double smc_calculate_probability(
//   vec& rankings,
//   uvec& unranked_items,
//   const vec& rho,
//   const double& alpha,
//   const bool& forward,
//   const std::string& metric = "footrule"
// ) {
//   double ranking_prob = 1;
//   for (uword jj = 0; jj < unranked_items.n_elem; ++jj) {
//
//
//     // next we get the sample probabilites for selecting a particular rank for
//     // an item based on the current alpha and the rho rank for that item
// //
// //     vec sample_prob_list = get_sample_probabilities(
// //       rho(unranked_items(jj)), alpha, rankings(unranked_items), rho.n_elem, metric
// //     );
//
//     // vec sampled_rank;
//     // if(forward){
//     //   sampled_rank = sample(rankings(unranked_items), 1, false, sample_prob_list);
//     // }
//     //
//     // double sampled_prob = as_scalar(sample_prob_list(find(rankings == sampled_rank)));
//     // ranking_prob *= sampled_prob;
//     //
//     // unranked_items = unranked_items(find(rankings != sampled_rank));
//   }
//
//   return ranking_prob;
// }
//
// double calculate_backward_probability(
//     arma::uvec item_ordering,
//     arma::vec partial_ranking,
//     arma::vec current_ranking,
//     arma::vec remaining_set,
//     const arma::vec rho,
//     const double alpha,
//     const int n_items,
//     const std::string metric = "footrule"
// ) {
//   // given an old and new item ordering, sample ranking with new ordering and
//   // calc the forward and backward prob
//
//   // show the augmented parts of the current ranking
//   // item ordering is the order of which items are assigned ranks in a specified
//   // order
//   const uword& num_items_unranked = item_ordering.n_elem;
//
//   // initialise prob of creating augmented ranking
//   double backward_auxiliary_ranking_probability = 1.0;
//
//   if (num_items_unranked != 1) {
//     // show the augmented parts of the current ranking
//     current_ranking = current_ranking.elem(item_ordering);
//
//     /* ====================================================== */
//     /* LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY     */
//     /* ====================================================== */
//     // given the old and new item ordering and the list of missing rank, determine
//     // the sample probs for each iteration
//     // backward_auxiliary_ranking_probability =
//     // smc_calculate_probability(
//     //   backward_auxiliary_ranking_probability, remaining_set, item_ordering,
//     //   current_ranking, rho, num_items_unranked, alpha, n_items, false, metric
//     // );
//   }
//   return 1;
// }
//
//
// Rcpp::List calculate_forward_probability(
//     vec& rankings,
//     uvec& unranked_items,
//     const double& alpha, const vec& rho,
//     const std::string& metric
// ) {
//
//   // prob of creating augmented ranking
//   double forward_auxiliary_ranking_probability = 1.0;
//
//
//   /* ====================================================== */
//   /* LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY     */
//   /* ====================================================== */
//   // given the old and new item ordering and the list of missing rank,
//   // determine the sample probs for each iteration
//
//   forward_auxiliary_ranking_probability = smc_calculate_probability(
//     rankings, unranked_items, rho, alpha, true, metric
//   );
//
//   // fit the augmented ranking within the partial rankings with NAs
//   // const vec& ar = conv_to<vec>::from(auxiliary_ranking);
//   // partial_ranking.elem(item_ordering) = ar; // ranks for items
//
//   return Rcpp::List::create(
//     Rcpp::Named("aug_ranking") = rankings,
//     Rcpp::Named("forward_prob") = forward_auxiliary_ranking_probability
//   );
// }
