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
  return distfun->matdist(rankings, rho);
}

double CayleyDistance::d(const arma::vec& r1, const arma::vec& r2) {
  double distance{};
  arma::vec tmp2 = r1;

  for(size_t i{}; i < r1.n_elem; ++i){
    if(tmp2(i) != r2(i)) {
      distance += 1;
      double tmp1 = tmp2(i);
      tmp2(i) = r2(i);
      arma::uvec inds = find(tmp2 == r2(i));
      tmp2.elem(inds).fill(tmp1);
    }
  }
  return distance;
}

arma::vec Distance::matdist(const arma::mat& r1, const arma::vec& r2) {
  arma::vec result(r1.n_cols);
  for(size_t i{}; i < r1.n_cols; i++) {
    const arma::vec& v1 = r1.col(i);
    result(i) = d(v1, r2);
  }
  return result;
}

arma::vec Distance::scalardist(const arma::vec& r1, const double r2) {
  arma::vec v2 = arma::ones(r1.size()) * r2;
  arma::vec result = arma::zeros(r1.size());
  for(size_t i{}; i < r1.n_elem; i++) {
    result(i) = d(arma::vec{r1(i)}, arma::vec{v2(i)});
  }
  return result;
}

void Distance::update_leap_and_shift_indices(arma::uvec& indices, int n_items) {
  return;
};

double CayleyDistance::d(
    const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1, r2);
}

void CayleyDistance::update_leap_and_shift_indices(
    arma::uvec& indices, int n_items) {
  indices = arma::regspace<arma::uvec>(0, n_items - 1);
}

double FootruleDistance::d(const arma::vec& r1, const arma::vec& r2) {
  return arma::norm(r1 - r2, 1);
}

double FootruleDistance::d(
    const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1(inds), r2(inds));
}

double HammingDistance::d(const arma::vec& r1, const arma::vec& r2) {
  return arma::sum(r1 != r2);
}

double HammingDistance::d(
    const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1(inds), r2(inds));
}

double KendallDistance::d(const arma::vec& r1, const arma::vec& r2) {
  double distance{};
  for(size_t i{}; i < r1.n_elem; ++i){
    for(size_t j{}; j < i; ++j){
      if (((r1(j) > r1(i)) && (r2(j) < r2(i)) ) ||
          ((r1(j) < r1(i)) && (r2(j) > r2(i)))) {
        distance += 1;
      }
    }
  }
  return distance;
}

double KendallDistance::d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1, r2);
}

double SpearmanDistance::d(const arma::vec& r1, const arma::vec& r2) {
  return std::pow(arma::norm(r1 - r2, 2), 2);
}

double SpearmanDistance::d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1(inds), r2(inds));
}

double UlamDistance::d(const arma::vec& r1, const arma::vec& r2) {
  arma::ivec a = arma::conv_to<arma::ivec>::from(r1) - 1;
  arma::ivec b = arma::conv_to<arma::ivec>::from(r2) - 1;
  auto distance = perm0_distance ( a, b );
  return static_cast<double>(distance);
}

double UlamDistance::d(const arma::vec& r1, const arma::vec& r2, const arma::uvec& inds) {
  return d(r1, r2);
}

void UlamDistance::update_leap_and_shift_indices(arma::uvec& indices, int n_items) {
  indices = arma::regspace<arma::uvec>(0, n_items - 1);
}
