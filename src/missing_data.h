#ifndef MISSING_H
#define MISSING_H

#include <RcppArmadillo.h>

void define_missingness(arma::mat& missing_indicator, arma::vec& assessor_missing,
                        const arma::mat& R,
                        const int& n_items, const int& n_assessors);
#endif
