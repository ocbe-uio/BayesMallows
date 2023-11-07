#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

arma::vec normalize_weights(const arma::vec& log_inc_wgt);
arma::vec initialize_alpha(const int& n_particles, const Rcpp::Nullable<arma::vec>& alpha_init = R_NilValue);
double get_exponent_sum(double, arma::vec, int, arma::mat, std::string);

arma::vec get_sample_probabilities(const arma::vec rho_item_rank,
                                   const double alpha,
                                   const arma::vec remaining_set_ranks,
                                   const int n_items,
                                   const std::string metric);
Rcpp::List calculate_forward_probability(const arma::vec&, const arma::uvec&, const arma::uvec&, const double&, const arma::vec&, const std::string&);
double calculate_backward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, arma::vec, double, int, std::string);
Rcpp::List correction_kernel(arma::vec, arma::vec, int);
Rcpp::List correction_kernel_pseudo(arma::vec, arma::vec, arma::vec, double, int, std::string);
arma::vec metropolis_hastings_aug_ranking(const double&, const arma::vec&, const int&, const arma::vec&, const arma::vec&, const bool&, const arma::uvec&, const std::string&);
#endif
