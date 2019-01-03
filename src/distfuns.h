#ifndef DISTFUNS_H
#define DISTFUNS_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "subset.h"


double get_rank_distance(arma::vec, arma::vec, std::string);
double rank_dist_matrix(const arma::mat&, const arma::vec&, std::string);

arma::vec update_distance_matrix(const arma::mat& rankings, const arma::vec& rho_cluster,
                                 const std::string& metric);

#endif
