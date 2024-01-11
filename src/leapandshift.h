#pragma once
#include <RcppArmadillo.h>
#include <memory>
#include "distances.h"
void leap_and_shift(arma::vec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::vec& rho, int leap_size,
                    const std::unique_ptr<Distance>& distfun);
void shift_step(arma::vec& rho_proposal, const arma::vec& rho,
                const int& u, arma::uvec& indices);
