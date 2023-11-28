#ifndef MISC_H
#define MISC_H

#include "RcppArmadillo.h"

double rtruncbeta(int shape1, int shape2, double trunc = 1);
arma::vec normalize_weights(const arma::vec& log_inc_wgt);

#endif
