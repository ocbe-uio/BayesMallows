#include <RcppArmadillo.h>
#include "misc.h"
#include "distfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]

void define_missingness(arma::mat& missing_indicator, arma::vec& assessor_missing,
                        const arma::mat& R,
                        const int& n_items, const int& n_assessors){
  for(int i = 0; i < n_assessors; ++i){
    for(int j = 0; j < n_items; ++j){
      if(!arma::is_finite(R(j, i))){
        missing_indicator(j, i) = 1;
        ++assessor_missing(i);
      }
    }
  }

}


arma::vec propose_augmentation(const arma::vec& ranks, const arma::vec& indicator,
                               const int& n_items){

  // Defining a number of helper variables
  arma::uvec ranked_inds = arma::find(indicator == 0);
  arma::uvec taken_ranks = arma::conv_to<arma::uvec>::from(ranks).elem(ranked_inds);
  arma::uvec all_ranks = arma::regspace<arma::uvec>(1, 1, n_items);
  arma::uvec all_inds = arma::regspace<arma::uvec>(0, 1, n_items - 1);

  // Find the available ranks
  arma::uvec available_ranks = std_setdiff(all_ranks, taken_ranks);
  arma::uvec available_inds = std_setdiff(all_inds, ranked_inds);

  // Adding randomness with shuffle
  arma::vec proposal = ranks;
  proposal(available_inds) = arma::shuffle(arma::conv_to<arma::vec>::from(available_ranks));

  return proposal;


}


void initialize_missing_ranks(arma::mat& R, const arma::mat& missing_indicator,
                              const arma::vec& assessor_missing,
                              const int& n_items, const int& n_assessors) {

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0) {
      continue;
    } else {
      R.col(i) = propose_augmentation(R.col(i), missing_indicator.col(i), n_items);
    }
  }
}

void update_missing_ranks(arma::mat& R, arma::mat& aug_acceptance,
                          const arma::mat& missing_indicator,
                          const arma::vec& assessor_missing,
                          const int& n_items, const int& n_assessors,
                          const double& alpha, const arma::vec& rho,
                          const std::string& metric, const int& t,
                          int& aug_diag_index, const int& aug_diag_thinning){
  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0){
      ++aug_acceptance(i, aug_diag_index);
    } else {
      arma::vec proposal = propose_augmentation(R.col(i), missing_indicator.col(i),
                                                n_items);

      // Draw a uniform random number
      double u = log(arma::randu<double>());

      double ratio = -alpha / n_items *
        (get_rank_distance(proposal, rho, metric) -
        get_rank_distance(R.col(i), rho, metric));

      if(ratio > u){
        R.col(i) = proposal;
        ++aug_acceptance(i, aug_diag_index);
      }
    }
  }

  // If appropriate, increment the augmentation diagnostic index
  if(t % aug_diag_thinning == 0){
    ++aug_diag_index;
  }

}
