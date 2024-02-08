#pragma once
#include <RcppArmadillo.h>

struct Resampler {
  Resampler() {};
  virtual ~Resampler() = default;
  virtual Rcpp::IntegerVector resample(arma::vec probs) = 0;
};

struct Multinomial : Resampler {
  Rcpp::IntegerVector resample(arma::vec probs) override;
};

struct Residual : Resampler {
  Rcpp::IntegerVector resample(arma::vec probs) override;
};

struct Stratified : Resampler {
  Rcpp::IntegerVector resample(arma::vec probs) override;
};

struct Systematic : Resampler {
  Rcpp::IntegerVector resample(arma::vec probs) override;
};
