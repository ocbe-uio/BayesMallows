#ifndef PARTITIONFUNS_H
#define PARTITIONFUNS_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "distfuns.h"

double get_partition_function(int, double, Rcpp::Nullable<arma::vec>,
                              Rcpp::Nullable<arma::vec>, std::string);

arma::vec asymptotic_partition_function(arma::vec alpha_grid, std::string metric,
                                        int K, int n_iterations, int n_items);

#endif
