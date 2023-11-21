#ifndef PAIRWISE_H
#define PAIRWISE_H

arma::vec propose_pairwise_augmentation(const arma::vec& ranking, const Rcpp::List& assessor_constraints);
arma::vec propose_swap(const arma::vec& ranking, const Rcpp::List& assessor_constraints,
                 int& g_diff, const int& swap_leap);

#endif
