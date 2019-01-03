#ifndef PARTITIONFUNS_H
#define PARTITIONFUNS_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "distances.h"

double get_partition_function(int, double, Rcpp::Nullable<arma::vec>,
                              Rcpp::Nullable<arma::vec>, std::string);

#endif
