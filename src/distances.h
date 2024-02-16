#pragma once
#include <memory>
#include <RcppArmadillo.h>

int perm0_distance (
    const arma::ivec& a,
    const arma::ivec& b);

struct Distance {
  Distance() {};
  virtual ~Distance() = default;
  virtual double d(const arma::vec& r1, const arma::vec& r2) = 0;
  virtual double d(const arma::vec& r1, const arma::vec& r2,
                   const arma::uvec& inds) = 0;
  arma::vec matdist(const arma::mat& r1, const arma::vec& r2);
  arma::vec scalardist(const arma::vec& r1, const double r2);
  virtual void update_leap_and_shift_indices(arma::uvec& indices, int n_items);
};

std::unique_ptr<Distance> choose_distance_function(std::string metric);

struct CayleyDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override;
  double d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) override;
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override;
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
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override;
};
