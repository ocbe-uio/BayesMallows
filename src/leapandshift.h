#ifndef LEAPANDSHIFT_H
#define LEAPANDSHIFT_H

#include "RcppArmadillo.h"

void leap_and_shift(arma::vec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::vec& rho, int leap_size, bool reduce_indices);

void shift_step(arma::vec& rho_proposal, const arma::vec& rho,
                const int& u, double& delta_r, arma::uvec& indices);

#endif
