#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

arma::vec normalize_weights(const arma::vec& log_inc_wgt);
arma::vec initialize_alpha(const int& N);
double get_exponent_sum(double, arma::vec, int, arma::mat, std::string);
arma::vec metropolis_hastings_rho(double, int, arma::mat, std::string, arma::vec, int);
double metropolis_hastings_alpha(double, int, arma::mat, std::string, arma::vec, const Rcpp::Nullable<arma::vec>, double, double, double);
arma::vec get_sample_probabilities(arma::vec, double, arma::vec, std::string, int);
Rcpp::List calculate_forward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, double, int, std::string);
double calculate_backward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, arma::vec, double, int, std::string);
Rcpp::List correction_kernel(arma::vec, arma::vec, int);
Rcpp::List correction_kernel_pseudo(arma::vec, arma::vec, arma::vec, double, int, std::string);
arma::vec metropolis_hastings_aug_ranking(const double&, const arma::vec&, const int&, const arma::vec&, const arma::vec&, const std::string&, const bool&);

#endif
