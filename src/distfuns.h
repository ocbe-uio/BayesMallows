#ifndef DISTFUNS_H
#define DISTFUNS_H

#include "RcppArmadillo.h"
#include "misc.h"
arma::vec get_summation_distances(int n, arma::vec cardinalities,
                                  std::string metric = "footrule");

double get_rank_distance(arma::vec, arma::vec, std::string);
double rank_dist_matrix(const arma::mat&, const arma::vec&, std::string);

#endif
