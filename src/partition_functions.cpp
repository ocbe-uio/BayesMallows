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
