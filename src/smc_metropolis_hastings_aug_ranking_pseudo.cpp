#include "RcppArmadillo.h"
#include "smc.h"
#include "misc.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Augmented Ranking (pseudolikelihood)
//' @description Function to perform Metropolis-Hastings for new augmented ranking using the pseudolikelihood augmentation approach
//'
//' @param alpha Numeric value of the scale parameter
//' @param rho Numeric vector specifying the consensus ranking
//' @param n_items Integer is the number of items in a ranking
//' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obtained from calculate_forward_probability function
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return = proposed augmented ranking or current ranking A ranking sequence vector representing proposed augmented ranking for next
//'         iteration of MCMC chain
//' @export
// [[Rcpp::export]]

arma::vec metropolis_hastings_aug_ranking_pseudo(
	double alpha,
	arma::vec rho,
	int n_items,
	arma::vec partial_ranking,
	arma::vec current_ranking,
	std::string metric
) {
  // augment incomplete ranks to initialise
  arma::vec ranks;
  Rcpp::IntegerVector tmp = Rcpp::seq(1, n_items);
  ranks = Rcpp::as<arma::vec>(tmp);

  // find items missing from original observed ranking
  arma::uvec unranked_items = find_nonfinite(partial_ranking);

  // find unallocated ranks from original observed ranking
  // Rcpp::NumericVector p_rank_Rcpp, c_rank_Rcpp;
  // p_rank_Rcpp = partial_ranking;
  // c_rank_Rcpp = current_ranking;
  arma::vec remaining_set = arma_setdiff_vec(current_ranking, partial_ranking);

  // if the observed and augmented ranking are exactly the same then break
  bool condition_1 = arma::approx_equal(\
    partial_ranking, current_ranking, "absdiff", 0.1\
  );
  bool condition_2 = remaining_set.n_elem == 1;
  if (condition_1 | condition_2) {
    return(current_ranking);
  } else {
    // randomly permute the unranked items to give the order in which they will
    // be allocated
    Rcpp::IntegerVector unranked_items_Rcpp;
    unranked_items_Rcpp = arma::conv_to<arma::ivec>::from(unranked_items);
    arma::uvec item_ordering = Rcpp::as<arma::uvec>(\
      Rcpp::sample(unranked_items_Rcpp, unranked_items_Rcpp.length())\
    );
    item_ordering = item_ordering + 1;

    // Calculate probabilities
    Rcpp::List proposal = calculate_forward_probability(\
      item_ordering, partial_ranking, remaining_set, rho, alpha, n_items,\
      metric\
    );
    arma::vec proposed_augmented_ranking = proposal["aug_ranking"];
    double forward_prob = proposal["forward_prob"];

    double backward_prob = calculate_backward_probability(\
      item_ordering, partial_ranking, current_ranking, remaining_set, rho,\
      alpha, n_items, metric\
    );

     /* MH TIME ------------------------------------------------------------- */
    // Calculate the log posterior of the current and proposed rankings.
    // NB the current can usually be stored to save recalculating it, but we're not caring about that yet
    double curr_logpost = get_mallows_loglik(\
      alpha, rho, n_items, current_ranking.t(), metric\
    );
    double prop_logpost = get_mallows_loglik(\
      alpha, rho, n_items, proposed_augmented_ranking.t(), metric\
    );

    double log_acceptance_prob = prop_logpost - curr_logpost - std::log(forward_prob) + std::log(backward_prob);

    if (std::log(Rcpp::as<double>(Rcpp::runif(1, 0, 1))) < log_acceptance_prob) {
      return(proposed_augmented_ranking);
    } else {
      return(current_ranking);
    }

  }
}
