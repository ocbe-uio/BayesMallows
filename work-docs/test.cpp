#include <memory>
#include "../src/partition_functions.h"
using namespace Rcpp;

// [[Rcpp::export]]
double logz(double alpha, int n_items, std::string metric,
            arma::vec distances, arma::vec cardinalities) {
  auto pfun = choose_partition_function(n_items, metric, distances, cardinalities);
  return pfun->logz(alpha);
}

/*** R
library(BayesMallows)
n_items <- 5
metric <- "ulam"
alpha <- 2
card <- get_cardinalities(n_items, metric)
logz(2, n_items, metric, card$distance, card$value)
*/
