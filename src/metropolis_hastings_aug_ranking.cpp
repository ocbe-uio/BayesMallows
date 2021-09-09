#include "RcppArmadillo.h"
#include "smc.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Augmented Ranking (pseudolikelihood)
//' @description Function to perform Metropolis-Hastings for new augmented ranking using the pseudolikelihood augmentation approach
//'
//' @param alpha Numeric value og the scale parameter
//' @param rho Numeric vector specifying the consensus ranking
//' @param n_items Integer is the number of items in a ranking
//' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obatined from calculate_forward_probability function
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return = proposed augmented ranking or current ranking A ranking sequence vector representing proposed augmented ranking for next
//'         iteration of MCMC chain
//' @export
// [[Rcpp::export]]

arma::vec metropolis_hastings_aug_ranking_CPP(
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
  arma::vec remaining_set = current_ranking.elem(find_nonfinite(partial_ranking));

  // if the observed and augmented ranking are exactly the same then break
  bool condition_1 = arma::approx_equal(\
    partial_ranking, current_ranking, "absdiff", 0.1\
  );
  bool condition_2 = remaining_set.n_elem == 1;
  if (condition_1 | condition_2) {
    return(current_ranking);
  } else {


    // generate random order for remaining_set
    remaining_set = std::move(arma::shuffle(remaining_set));

    // Subset by element position and set equal to the now permuted remaining set
    partial_ranking.elem(unranked_items) = remaining_set;

    // set the augmented partial ranking as the proposed augmented ranking
    arma::vec proposed_ranking = partial_ranking;


     /* MH TIME ------------------------------------------------------------- */
    // Calculate the log posterior of the current and proposed rankings.
    // NB the current can usually be stored to save recalculating it, but we're not caring about that yet
    double curr_logpost = get_mallows_loglik(\
      alpha, rho, n_items, current_ranking.t(), metric\
    );
    double prop_logpost = get_mallows_loglik(\
      alpha, rho, n_items, proposed_ranking.t(), metric\
    );

    double log_acceptance_prob = prop_logpost - curr_logpost;

    if (std::log(Rcpp::as<double>(Rcpp::runif(1, 0, 1))) < log_acceptance_prob) {
      return(proposed_ranking);
    } else {
      return(current_ranking);
    }

  }
}