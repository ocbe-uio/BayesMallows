#ifndef IMPORTANCESAMPLING_H
#define IMPORTANCESAMPLING_H

#include "RcppArmadillo.h"
#include "distfuns.h"

arma::vec compute_importance_sampling_estimate(arma::vec, int, std::string, int);

#endif
