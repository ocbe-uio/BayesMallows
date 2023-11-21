#ifndef MIXTURES_H
#define MIXTURES_H

#include "RcppArmadillo.h"

void update_dist_mat(arma::mat& dist_mat, const arma::mat& rankings,
                     const arma::mat& rho_old, const std::string& metric,
                     const arma::vec& observation_frequency);


#endif
