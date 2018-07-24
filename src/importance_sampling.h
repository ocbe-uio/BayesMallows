#ifndef IMPORTANCESAMPLING_H
#define IMPORTANCESAMPLING_H

#include "RcppArmadillo.h"
#include "distfuns.h"

arma::vec get_is_estimate(arma::vec, int, std::string, std::string, int);

#endif
