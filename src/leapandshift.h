#pragma once
#include <RcppArmadillo.h>
#include <memory>
#include "distances.h"

struct LeapShiftObject{
  LeapShiftObject(const arma::vec& rho_proposal,
                  const arma::uvec& indices = arma::uvec{},
                  double prob_backward = 1, double prob_forward = 1) :
  rho_proposal { rho_proposal },
  indices { indices },
  prob_backward { prob_backward },
  prob_forward { prob_forward }{};

  ~LeapShiftObject() = default;

  void shift_step(const arma::vec& rho, int u);

  arma::vec rho_proposal;
  arma::uvec indices;
  double prob_backward;
  double prob_forward;
};

LeapShiftObject leap_and_shift(
    const arma::vec& rho, int leap_size,
    const std::unique_ptr<Distance>& distfun);
