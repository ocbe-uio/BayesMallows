#ifndef DISTANCES_H
#define DISTANCES_H

#include "RcppArmadillo.h" // needed because of rank_dist_vec function


double get_rank_distance(arma::uvec, arma::uvec, std::string);
double rank_dist_sum(const arma::umat&, const arma::uvec&, const std::string&, const arma::uvec&);

arma::vec rank_dist_vec(const arma::umat& rankings, const arma::uvec& rho,
                        const std::string& metric, const arma::uvec& obs_freq);

#endif
