#include "distances.h"

std::unique_ptr<Distance> choose_distance_function(std::string metric) {
  if(metric == "cayley") {
    return std::unique_ptr<Distance>{new CayleyDistance{}};
  } else if(metric == "footrule") {
    return std::unique_ptr<Distance>{new FootruleDistance{}};
  } else if(metric == "hamming") {
    return std::unique_ptr<Distance>{new HammingDistance{}};
  } else if(metric == "kendall") {
    return std::unique_ptr<Distance>{new KendallDistance{}};
  } else if(metric == "spearman") {
    return std::unique_ptr<Distance>{new SpearmanDistance{}};
  } else if(metric == "ulam") {
    return std::unique_ptr<Distance>{new UlamDistance{}};
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
