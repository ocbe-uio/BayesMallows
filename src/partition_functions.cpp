#include "partition_functions.h"

std::unique_ptr<PartitionFunction> choose_partition_function(
    int n_items, std::string metric,
    const Rcpp::Nullable<arma::mat>& pfun_values,
    const Rcpp::Nullable<arma::mat>& pfun_estimate) {
  if(metric == "cayley") {
    return std::make_unique<Cayley>(n_items);
  } else if(metric == "hamming") {
    return std::make_unique<Hamming>(n_items);
  } else if(metric == "kendall") {
    return std::make_unique<Kendall>(n_items);
  } else if(pfun_values.isNotNull()) {
    return std::make_unique<Cardinal>(
      n_items, Rcpp::as<arma::mat>(pfun_values));
  } else if(pfun_estimate.isNotNull()) {
    return std::make_unique<Estimated>(
      n_items, Rcpp::as<arma::mat>(pfun_estimate));
  } else {
    Rcpp::stop("Unknown metric.");
  }
}

// [[Rcpp::export]]
double get_expected_distance(
    double alpha, int n_items, std::string metric,
    const Rcpp::Nullable<arma::mat>& pfun_values) {
  auto pfun = choose_partition_function(
    n_items, metric, pfun_values, R_NilValue);
  return pfun->expected_distance(alpha);
}

// [[Rcpp::export]]
double get_partition_function(
    double alpha, int n_items, std::string metric,
    const Rcpp::Nullable<arma::mat>& pfun_values) {
  auto pfun = choose_partition_function(
    n_items, metric, pfun_values, R_NilValue);
  return pfun->logz(alpha);
}

Cayley::Cayley(int n_items) : n_items { n_items } {}

double Cayley::logz(double alpha) {
  double res{};
  for(int i{1}; i < n_items; ++i){
    res += std::log(1.0 + i * std::exp(- alpha / n_items));
  }
  return res;
}

double Cayley::expected_distance(double alpha) {
  arma::vec idx = arma::regspace<arma::vec>(1, n_items - 1);
  return arma::sum(idx / (idx + std::exp(alpha / n_items)));
}

Hamming::Hamming(int n_items) : n_items { n_items } {}

double Hamming::logz(double alpha) {
  double res{};
  for(int i{}; i < (n_items + 1); ++i){
    res += tgamma(n_items + 1.0) * std::exp(-alpha) *
      std::pow((std::exp(alpha / n_items) - 1.0), i) / tgamma(i + 1.0);
  }
  return std::log(res);
}

double Hamming::expected_distance(double alpha) {
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

Kendall::Kendall(int n_items) : n_items { n_items } {}

double Kendall::logz(double alpha) {
  double res{};
  for(int i{1}; i < (n_items + 1); ++i){
    res += std::log((1.0 - std::exp(-i * alpha / n_items)) /
      (1.0 - std::exp(-alpha / n_items)));
  }
  return res;
}

double Kendall::expected_distance(double alpha) {
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

Cardinal::Cardinal(
  int n_items,
  const arma::mat& pfun_values) :
  n_items { n_items },
  distances { pfun_values.col(0) },
  cardinalities { pfun_values.col(1) }{}

double Cardinal::logz(double alpha) {
  return std::log(arma::sum(cardinalities % exp(-alpha / n_items * distances)));
}

double Cardinal::expected_distance(double alpha) {
  return arma::sum(distances % cardinalities % exp(-alpha * distances / n_items))
  * std::exp(-logz(alpha));
}

Estimated::Estimated(
  int n_items,
  const arma::mat& pfun_estimate) :
  n_items { n_items },
  power { pfun_estimate.col(0) },
  coefficients { pfun_estimate.col(1) } {}

double Estimated::logz(double alpha) {
  return arma::sum(
    arma::pow(alpha + arma::zeros(coefficients.size()), power) % coefficients
  );
}

double Estimated::expected_distance(double alpha) {
  Rcpp::stop(
    "Expected distance not available with estimated partition function.");
}
