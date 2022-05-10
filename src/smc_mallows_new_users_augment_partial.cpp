#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void smc_mallows_new_users_augment_partial(
    arma::cube& aug_rankings,
    arma::vec& aug_prob,
    const arma::cube& rho_samples,
    const arma::mat& alpha_samples,
    const int& num_obs,
    const int& num_new_obs,
    const arma::mat& R_obs,
    const std::string& aug_method,
    const int& N,
    const std::string& metric,
    const int& tt,
    const int& n_items,
    const double& alpha,
    const bool& augment_alpha
){
  ivec ranks = regspace<ivec>(1, n_items);
  for (int ii{}; ii < N; ++ii) {
    for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
      vec partial_ranking = R_obs.row(jj).t();

      // find items missing from original observed ranking
      const uvec& unranked_items = find_nonfinite(partial_ranking);

      // find ranks missing from ranking
      const vec& missing_ranks = setdiff_template(ranks, partial_ranking);

      // fill in missing ranks based on choice of augmentation method
      Rcpp::List proposal{};
      if (aug_method == "random") {
        // create new augmented ranking by sampling remaining ranks from set uniformly
        partial_ranking.elem(find_nonfinite(partial_ranking)) = shuffle(missing_ranks);

        aug_rankings(span(jj), span::all, span(ii)) = partial_ranking;
        aug_prob(ii) = divide_by_fact(aug_prob(ii), missing_ranks.n_elem);
      } else if ((aug_method == "pseudolikelihood") && ((metric == "footrule") || (metric == "spearman"))) {
        // randomly permute the unranked items to give the order in which they will be allocated
        uvec item_ordering = shuffle(unranked_items);
        const rowvec& rho_s = rho_samples(span(ii), span::all, span(tt + 1));
        proposal = calculate_forward_probability(
          item_ordering, partial_ranking, missing_ranks, rho_s.t(),
          augment_alpha ? alpha_samples(ii, tt + 1) : alpha, n_items, metric
        );

        const vec& a_rank = proposal["aug_ranking"];
        const double& f_prob = proposal["forward_prob"];
        aug_rankings(span(jj), span::all, span(ii)) = a_rank;
        aug_prob(ii) = aug_prob(ii) * f_prob;
      } else {
        Rcpp::stop("Combined choice of metric and aug_method is incompatible");
      }
    }
  }
}

