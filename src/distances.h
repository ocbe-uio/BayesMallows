#ifndef DISTANCES_H
#define DISTANCES_H

#include "RcppArmadillo.h"
#include "subset.h"


double get_rank_distance(arma::vec, arma::vec, std::string);
double rank_dist_sum(const arma::mat&, const arma::vec&, const std::string&);

arma::vec rank_dist_vec(const arma::mat& rankings, const arma::vec& rho, const std::string& metric);

#endif
