#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

arma::vec normalize_weights(const arma::vec& log_inc_wgt);
arma::vec initialize_alpha(const int& N);
double get_exponent_sum(double, arma::uvec, uint, arma::umat, std::string);
arma::uvec metropolis_hastings_rho(double, uint, arma::umat, arma::uvec, std::string, int);
double metropolis_hastings_alpha(double, uint, arma::umat, arma::uvec, const Rcpp::Nullable<arma::vec>, std::string, double, double, double);
arma::vec get_sample_probabilities(arma::uvec, double, arma::uvec, int, std::string);
Rcpp::List calculate_forward_probability(arma::uvec, arma::uvec, arma::uvec, arma::uvec, double, int, std::string);
double calculate_backward_probability(arma::uvec, arma::uvec, arma::uvec, arma::uvec, arma::uvec, double, int, std::string);
Rcpp::List correction_kernel(arma::uvec, arma::uvec, uint);
Rcpp::List correction_kernel_pseudo(arma::uvec, arma::uvec, arma::uvec, double, int, std::string);
arma::uvec metropolis_hastings_aug_ranking(const double&, const arma::uvec&, const uint&, const arma::uvec&, const arma::uvec&, const bool&, const std::string&);
#endif
