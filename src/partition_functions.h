#pragma once

#include <RcppArmadillo.h>

struct PartitionFunction {
  PartitionFunction() {};
  virtual ~PartitionFunction() = default;
  virtual double logz(double alpha) = 0;
  virtual double expected_distance(double alpha) = 0;
};

std::unique_ptr<PartitionFunction> choose_partition_function(
    int n_items, std::string metric, const arma::vec& distances,
    const arma::vec& cardinalities);

struct Cayley : PartitionFunction {
  Cayley(int n_items) : n_items { n_items } {}

  double logz(double alpha) override {
    double res{};
    for(int i{1}; i < n_items; ++i){
      res += std::log(1.0 + i * std::exp(- alpha / n_items));
    }
    return res;
  }
  double expected_distance(double alpha) override {
    arma::vec idx = arma::regspace<arma::vec>(1, n_items - 1);
    return arma::sum(idx / (idx + std::exp(alpha / n_items)));
  }
  const int n_items;
};

struct Hamming : PartitionFunction {
  Hamming(int n_items) : n_items { n_items } {}

  double logz(double alpha) override {
    double res{};
    for(int i{}; i < (n_items + 1); ++i){
      res += tgamma(n_items + 1.0) * std::exp(-alpha) *
        std::pow((std::exp(alpha / n_items) - 1.0), i) / tgamma(i + 1.0);
    }
    return std::log(res);
  }
  double expected_distance(double alpha) override {
    arma::vec idx1 = arma::regspace<arma::vec>(0, n_items - 1);
    arma::vec idx2 = arma::regspace<arma::vec>(0, n_items);

    return n_items - std::exp(alpha / n_items) *
      arma::sum(
        arma::pow(std::exp(alpha / n_items) - arma::ones(idx1.size()), idx1) /
        arma::tgamma(idx1 + 1)
      ) /
        arma::sum(
          arma::pow(std::exp(alpha / n_items) - arma::ones(idx2.size()), idx2) /
          arma::tgamma(idx2 + 1));
  }

  const int n_items;
};

struct Kendall : PartitionFunction {
  Kendall(int n_items) : n_items { n_items } {}

  double logz(double alpha) override {
    double res{};
    for(int i{1}; i < (n_items + 1); ++i){
      res += std::log((1.0 - std::exp(-i * alpha / n_items)) /
        (1.0 - std::exp(-alpha / n_items)));
    }
    return res;
  }
  double expected_distance(double alpha) override {
    arma::vec idx = arma::regspace<arma::vec>(1, n_items);
    if(alpha > 0) {
      return n_items * std::exp(-alpha / n_items) /
        (1 - std::exp(-alpha / n_items)) -
          arma::sum((idx % arma::exp(-idx * alpha / n_items)) /
            (1 - arma::exp(-idx * alpha / n_items)));
    } else if(alpha == 0) {
      return n_items * (n_items - 1) / 4;
    } else {
      Rcpp::stop("alpha must be non-negative.");
    }
  }

  const int n_items;
};

struct Cardinal : PartitionFunction {
  Cardinal(
    int n_items,
    const arma::vec& distances,
    const arma::vec& cardinalities) :
  n_items { n_items },
  distances { distances },
  cardinalities { cardinalities }{}

  double logz(double alpha) override {
    return std::log(arma::sum(cardinalities % exp(-alpha / n_items * distances)));
  }
  double expected_distance(double alpha) override {
    return arma::sum(distances % cardinalities % exp(-alpha * distances / n_items))
    * std::exp(-logz(alpha));
  }

  const int n_items;
  const arma::vec distances;
  const arma::vec cardinalities;
};

