#include "distances.h"
using namespace arma;

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
arma::vec get_rank_distance(
    arma::mat rankings, arma::vec rho, std::string metric) {
  auto distfun = choose_distance_function(metric);
  return distfun->matdist(rankings, rho);
}

vec Distance::matdist(const mat& r1, const vec& r2) {
  return matdist(r1, r2, regspace<uvec>(0, r2.n_elem - 1));
}

vec Distance::matdist(const mat& r1, const vec& r2,
                            const uvec& inds) {
  vec result(r1.n_cols);
  for(size_t i{}; i < r1.n_cols; i++) {
    const vec& v1 = r1.col(i);
    result(i) = d(v1, r2, inds);
  }
  return result;
}

vec Distance::scalardist(const vec& r1, const double r2) {
  vec v2 = ones(r1.size()) * r2;
  vec result = zeros(r1.size());
  for(size_t i{}; i < r1.n_elem; i++) {
    result(i) = d(vec{r1(i)}, vec{v2(i)});
  }
  return result;
}

double CayleyDistance::d(const vec& r1, const vec& r2) {
  double distance{};
  vec tmp2 = r1;

  for(size_t i{}; i < r1.n_elem; ++i){
    if(tmp2(i) != r2(i)) {
      distance += 1;
      double tmp1 = tmp2(i);
      tmp2(i) = r2(i);
      uvec inds = find(tmp2 == r2(i));
      tmp2.elem(inds).fill(tmp1);
    }
  }
  return distance;
}

double CayleyDistance::d(
    const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1, r2);
}

double FootruleDistance::d(const vec& r1, const vec& r2) {
  return norm(r1 - r2, 1);
}

double FootruleDistance::d(
    const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1(inds), r2(inds));
}

double HammingDistance::d(const vec& r1, const vec& r2) {
  return sum(r1 != r2);
}

double HammingDistance::d(
    const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1(inds), r2(inds));
}

double KendallDistance::d(const vec& r1, const vec& r2) {
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

double KendallDistance::d(const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1, r2);
}

double SpearmanDistance::d(const vec& r1, const vec& r2) {
  return std::pow(norm(r1 - r2, 2), 2);
}

double SpearmanDistance::d(const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1(inds), r2(inds));
}

// Rewritten from https://www.geeksforgeeks.org/c-program-for-longest-increasing-subsequence/
int longest_increasing_subsequence(const vec& permutation) {
  int n = permutation.n_elem;
  vec lis(n, fill::ones);

  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      if (permutation(i) > permutation(j) && lis(i) < lis(j) + 1) {
        lis(i) = lis(j) + 1;
      }
    }
  }

  int max = 0;
  for (int i = 0; i < n; i++) {
    if (max < lis(i)) {
      max = lis(i);
    }
  }

  return max;
}


double UlamDistance::d(const vec& r1, const vec& r2) {
  uvec x = sort_index(r2);
  return r1.size() - longest_increasing_subsequence(r1(x));
}

double UlamDistance::d(const vec& r1, const vec& r2, const uvec& inds) {
  return d(r1, r2);
}
