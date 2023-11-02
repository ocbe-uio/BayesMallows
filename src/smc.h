#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

arma::vec normalize_weights(const arma::vec& log_inc_wgt);
arma::vec initialize_alpha(const int& n_particles, const Rcpp::Nullable<arma::vec>& alpha_init = R_NilValue);
double get_exponent_sum(double, arma::vec, int, arma::mat, std::string);
arma::vec metropolis_hastings_rho(double, int, arma::mat, arma::vec, std::string, int);
double metropolis_hastings_alpha(double, int, arma::mat, arma::vec, const Rcpp::Nullable<arma::vec>, const Rcpp::Nullable<arma::vec>, std::string, double, double);
arma::vec get_sample_probabilities(arma::vec, double, arma::vec, int, std::string);
Rcpp::List calculate_forward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, double, int, std::string);
double calculate_backward_probability(arma::uvec, arma::vec, arma::vec, arma::vec, arma::vec, double, int, std::string);
Rcpp::List correction_kernel(arma::vec, arma::vec, int);
Rcpp::List correction_kernel_pseudo(arma::vec, arma::vec, arma::vec, double, int, std::string);
arma::vec metropolis_hastings_aug_ranking(const double&, const arma::vec&, const int&, const arma::vec&, const arma::vec&, const bool&, const std::string&);
#endif
