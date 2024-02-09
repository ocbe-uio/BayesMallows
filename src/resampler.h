#pragma once
#include <RcppArmadillo.h>

struct Resampler {
  Resampler() {};
  virtual ~Resampler() = default;
  virtual arma::ivec resample(arma::vec probs) = 0;
};

struct Multinomial : Resampler {
  arma::ivec resample(arma::vec probs) override;
};

struct Residual : Resampler {
  arma::ivec resample(arma::vec probs) override;
};

struct Stratified : Resampler {
  arma::ivec resample(arma::vec probs) override;
};

struct Systematic : Resampler {
  arma::ivec resample(arma::vec probs) override;
};

std::unique_ptr<Resampler> choose_resampler(std::string resampler);
