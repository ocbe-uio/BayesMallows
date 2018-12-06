#ifndef DISTFUNS_H
#define DISTFUNS_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "subset.h"


double get_rank_distance(arma::vec, arma::vec, std::string);
double rank_dist_matrix(const arma::mat&, const arma::vec&, std::string);

void update_distance_matrix(arma::mat& dist_mat, const arma::mat& rankings, const arma::vec& rho_cluster,
                            const int& n_assessors, const int& cluster_index,
                            const std::string metric);

#endif
