#include <RcppArmadillo.h>
#include "smc.h"
#include "misc.h"
#include "setdiff.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
arma::vec metropolis_hastings_aug_ranking_both(
    double alpha,
    arma::vec rho,
    int n_items,
    arma::vec partial_ranking,
    arma::vec current_ranking,
    std::string metric,
    bool pseudo
) {
  double forward_backward_prob{};
  vec proposed_ranking;
  // augment incomplete ranks to initialise
  vec ranks = regspace(1, n_items);

  // find items missing from original observed ranking
  const uvec unranked_items = find_nonfinite(partial_ranking);

  // find unallocated ranks from original observed ranking
  const vec remaining_set = setdiff_template(current_ranking, partial_ranking);

  // if the observed and augmented ranking are exactly the same then break
  const bool condition_1 = approx_equal(
    partial_ranking, current_ranking, "absdiff", 0.1
  );
  const bool condition_2 = remaining_set.n_elem == 1;
  if (condition_1 | condition_2) {
    return current_ranking;
  } else {
    if(pseudo){
      // randomly permute the unranked items to give the order in which they will
      // be allocated
      uvec item_ordering = shuffle(unranked_items);

      // Calculate probabilities
      const Rcpp::List proposal = calculate_forward_probability(
        item_ordering, partial_ranking, remaining_set, rho, alpha, n_items,
        metric
      );
      const vec& proposed_augmented_ranking = proposal["aug_ranking"];
      double forward_prob = proposal["forward_prob"];

      double backward_prob = calculate_backward_probability(
        item_ordering, partial_ranking, current_ranking, remaining_set, rho,
        alpha, n_items, metric
      );
      forward_backward_prob = - std::log(forward_prob) + std::log(backward_prob);
      proposed_ranking = proposed_augmented_ranking;

    } else {
      // Subset by element position and set equal to the now permuted remaining set
      partial_ranking.elem(unranked_items) = shuffle(remaining_set);

      // set the augmented partial ranking as the proposed augmented ranking
      proposed_ranking = partial_ranking;
    }
  }

  /* MH TIME ------------------------------------------------------------- */
  // Calculate the log posterior of the current and proposed rankings.
  // NB the current can usually be stored to save recalculating it, but we're not caring about that yet
  const double curr_logpost = get_exponent_sum(     \
    alpha, rho, n_items, current_ranking.t(), metric\
  );
  const double prop_logpost = get_exponent_sum(                \
    alpha, rho, n_items, proposed_ranking.t(), metric\
  );

  const double log_acceptance_prob = prop_logpost - curr_logpost + forward_backward_prob;

  if (std::log(Rcpp::as<double>(Rcpp::runif(1, 0, 1))) < log_acceptance_prob) {
    return(proposed_ranking);
  } else {
    return(current_ranking);
  }
}
