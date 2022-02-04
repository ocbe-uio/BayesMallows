#include "RcppArmadillo.h"
#include "misc.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Correction Kernel
//' @description Function to determine if the augmented ranking is compatible
//' with the new observed partial ranking. If it is not, the we create a new
//' augmentation using the random sampling approach and calculate the
//' augmentation probability.
//'
//' @param current_ranking A ranking sequence vector of the current augmented
//' ranking (no missing values)
//' @param observed_ranking A ranking sequence vector of the observed partial
//' ranking (no missing values) The original incomplete partial ranking
//' is in the rankings data set.
//' @param n_items Integer is the number of items in a ranking
//'
//' @return List containing the proposed 'corrected' augmented ranking
//' that is compatible with the new observed ranking for a user
// [[Rcpp::export]]
Rcpp::List correction_kernel(
  arma::vec observed_ranking, //R_obs
  arma::vec current_ranking,  // R_curr
  int n_items
) {
  // check if new information means 'mistakes' made with augmented rankings
  bool observed_equals_current = arma::approx_equal(\
    observed_ranking, current_ranking, "absdiff", 0.1\
  );
  double correction_prob = 1.0;
  arma::vec proposed_ranking;
  if (observed_equals_current) {
    proposed_ranking = current_ranking;
  } else {
    // resample from smaller pool of possible augmented rankings

    //  find elements missing from original observed ranking
    arma::vec remaining_set = arma_setdiff_vec(current_ranking, observed_ranking);

    // create new agumented ranking by sampling remaining ranks from set uniformly
    proposed_ranking = observed_ranking;

    arma::uvec unranked_items = find_nonfinite(proposed_ranking);
    if (remaining_set.n_elem == 1) {
      proposed_ranking.elem(unranked_items) = remaining_set;
    } else {
      // generate random order for remaining_set
      remaining_set = std::move(arma::shuffle(remaining_set));
      proposed_ranking.elem(unranked_items) = remaining_set;
    }

    Rcpp::NumericVector remaining_set_Rcpp;
    remaining_set_Rcpp = remaining_set;
    Rcpp::NumericVector remaining_set_Rcpp_elem;
    remaining_set_Rcpp_elem = remaining_set_Rcpp.length();
    Rcpp::NumericVector remaining_set_fact = Rcpp::factorial(remaining_set_Rcpp_elem);
    double remaining_set_fact_dbl = Rcpp::as<double>(remaining_set_fact);
    correction_prob = 1.0 / remaining_set_fact_dbl;
  }
  return Rcpp::List::create(
      Rcpp::Named("ranking") = proposed_ranking,
      Rcpp::Named("correction_prob") = correction_prob
  );
}
