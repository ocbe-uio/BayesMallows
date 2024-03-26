#pragma once
#include <memory>
#include <RcppArmadillo.h>

struct Distance {
  Distance() {};
  virtual ~Distance() = default;
  virtual double d(const arma::vec& r1, const arma::vec& r2) = 0;
  virtual double d(const arma::vec& r1, const arma::vec& r2,
                   const arma::uvec& inds) = 0;
  arma::vec matdist(const arma::mat& r1, const arma::vec& r2);
  arma::vec matdist(const arma::mat& r1, const arma::vec& r2,
                    const arma::uvec& inds);
  arma::vec scalardist(const arma::vec& r1, const double r2);
};

std::unique_ptr<Distance> choose_distance_function(std::string metric);

struct CayleyDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};

struct FootruleDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};

struct HammingDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};

struct KendallDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};

struct SpearmanDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};

struct UlamDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
};
