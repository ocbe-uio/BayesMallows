#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

double get_mallows_loglik(double, arma::vec, int, arma::mat, std::string);
Rcpp::List leap_and_shift_probs(arma::vec, int, int);
arma::vec metropolis_hastings_rho(double, int, arma::mat, std::string, arma::vec, int);
double metropolis_hastings_alpha(double, int, arma::mat, std::string, arma::vec, const Rcpp::Nullable<arma::vec>, double, double, double);
arma::vec get_sample_probabilities(arma::vec, double, arma::vec, std::string, int);
Rcpp::List calculate_forward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, double, int, std::string);
double calculate_backward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, arma::vec, double, int, std::string);

#endif
