#pragma once

#include <memory>
#include <RcppArmadillo.h>

struct PartitionFunction {
  PartitionFunction() {};
  virtual ~PartitionFunction() = default;
  virtual double logz(double alpha) = 0;
  virtual double expected_distance(double alpha) = 0;
};

std::unique_ptr<PartitionFunction> choose_partition_function(
    int n_items, std::string metric,
    const Rcpp::Nullable<arma::mat>& pfun_values,
    const Rcpp::Nullable<arma::mat>& pfun_estimate);

struct Cayley : PartitionFunction {
  Cayley(int n_items);
  double logz(double alpha) override;
  double expected_distance(double alpha) override;
  const int n_items;
};

struct Hamming : PartitionFunction {
  Hamming(int n_items);
  double logz(double alpha) override;
  double expected_distance(double alpha) override;
  const int n_items;
};

struct Kendall : PartitionFunction {
  Kendall(int n_items);
  double logz(double alpha) override;
  double expected_distance(double alpha) override;
  const int n_items;
};

struct Cardinal : PartitionFunction {
  Cardinal(int n_items, const arma::mat& pfun_values);
  double logz(double alpha) override;
  double expected_distance(double alpha) override;
  const int n_items;
  const arma::vec distances;
  const arma::vec cardinalities;
};

struct Estimated : PartitionFunction {
  Estimated(int n_items, const arma::mat& pfun_estimate);
  double logz(double alpha) override;
  double expected_distance(double alpha) override;
  const int n_items;
  const arma::vec power;
  const arma::vec coefficients;
};

