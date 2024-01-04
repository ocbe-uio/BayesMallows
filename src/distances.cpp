#include "distances.h"

std::unique_ptr<Distance> choose_distance_function(std::string metric) {
  if(metric == "cayley") {
    return std::make_unique<CayleyDistance>();
  } else if(metric == "footrule") {
    return std::make_unique<FootruleDistance>();
  } else if(metric == "hamming") {
    return std::make_unique<HammingDistance>();
  } else if(metric == "kendall") {
    return std::make_unique<KendallDistance>();
  } else if(metric == "spearman") {
    return std::make_unique<SpearmanDistance>();
  } else if(metric == "ulam") {
    return std::make_unique<UlamDistance>();
  } else {
    Rcpp::stop("Unknown metric.");
  }
}

// [[Rcpp::export]]
arma::vec get_rank_distance(arma::mat rankings, arma::vec rho,
                           std::string metric) {
  auto distfun = choose_distance_function(metric);
  return distfun->d(rankings, rho);
}
