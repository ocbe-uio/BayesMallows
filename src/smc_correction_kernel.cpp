// #include <RcppArmadillo.h>
// #include "misc.h"
// #include "setdiff.h"
//
// // [[Rcpp::depends(RcppArmadillo)]]
//
// Rcpp::List correction_kernel(
//   const arma::vec observed_ranking,
//   const arma::vec current_ranking,
//   const int n_items
// ) {
//   // check if new information means 'mistakes' made with augmented rankings
//   const bool observed_equals_current = approx_equal(\
//     observed_ranking, current_ranking, "absdiff", 0.1\
//   );
//   double correction_prob = 1.0;
//   arma::vec proposed_ranking;
//   if (observed_equals_current) {
//     proposed_ranking = current_ranking;
//   } else {
//     // resample from smaller pool of possible augmented rankings
//
//     //  find elements missing from original observed ranking
//     arma::vec remaining_set = setdiff_template(current_ranking, observed_ranking);
//
//     // create new agumented ranking by sampling remaining ranks from set uniformly
//     proposed_ranking = observed_ranking;
//
//     const arma::uvec unranked_items = find_nonfinite(proposed_ranking);
//     if (remaining_set.n_elem == 1) {
//       proposed_ranking.elem(unranked_items) = remaining_set;
//     } else {
//       // generate random order for remaining_set
//       remaining_set = std::move(shuffle(remaining_set));
//       proposed_ranking.elem(unranked_items) = remaining_set;
//     }
//     correction_prob = divide_by_fact(correction_prob, remaining_set.n_elem);
//   }
//   return Rcpp::List::create(
//       Rcpp::Named("ranking") = proposed_ranking,
//       Rcpp::Named("correction_prob") = correction_prob
//   );
// }
