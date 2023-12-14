#include <memory>
#include "partition_functions.h"

std::unique_ptr<PartitionFunction> choose_partition_function(
    int n_items, std::string metric, const arma::vec& distances,
    const arma::vec& cardinalities) {
  if(metric == "cayley") {
    return std::unique_ptr<PartitionFunction>{new Cayley{n_items}};
  } else if(metric == "hamming") {
    return std::unique_ptr<PartitionFunction>{new Hamming{n_items}};
  } else if(metric == "kendall") {
    return std::unique_ptr<PartitionFunction>{new Kendall{n_items}};
  } else if(metric == "footrule" || metric == "spearman" || metric == "ulam") {
    return std::unique_ptr<PartitionFunction>{
      new Cardinal{n_items, distances, cardinalities}};
  } else {
    Rcpp::stop("Unknown metric.");
  }
}

// [[Rcpp::export]]
double get_expected_distance(
    double alpha, int n_items, std::string metric,
    const arma::vec& distances,
    const arma::vec& cardinalities) {
  auto pfun = choose_partition_function(
    n_items, metric, distances, cardinalities);
  return pfun->expected_distance(alpha);
}

// [[Rcpp::export]]
double get_partition_function(
    double alpha, int n_items, std::string metric,
    const arma::vec& distances,
    const arma::vec& cardinalities) {
  auto pfun = choose_partition_function(
    n_items, metric, distances, cardinalities);
  return pfun->logz(alpha);
}
