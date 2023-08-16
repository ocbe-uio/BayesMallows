#ifndef LEAPANDSHIFT_H
#define LEAPANDSHIFT_H

#include <RcppArmadillo.h>

void leap_and_shift(arma::uvec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::uvec& rho, int leap_size, bool reduce_indices);

void shift_step(arma::uvec& rho_proposal, const arma::uvec& rho,
                const int& u, arma::uvec& indices);

#endif
