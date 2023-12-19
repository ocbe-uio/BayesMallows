#pragma once
#include <memory>
#include <RcppArmadillo.h>

int perm0_distance (
    const arma::ivec& a,
    const arma::ivec& b);

struct Distance {
  Distance() {};
  virtual ~Distance() = default;
  virtual double d(const arma::vec& r1, const arma::vec& r2) = 0;
  arma::vec d(const arma::mat& r1, const arma::vec& r2) {
    arma::vec result(r1.n_cols);
    for(size_t i{}; i < r1.n_cols; i++) {
      const arma::vec& v1 = r1.col(i);
      result(i) = d(v1, r2);
    }
    return result;
  }
  arma::vec d(const arma::vec& r1, const double r2) {
    arma::vec v2 = arma::ones(r1.size()) * r2;
    arma::vec result = arma::zeros(r1.size());
    for(size_t i{}; i < r1.n_elem; i++) {
      result(i) = d(arma::vec{r1(i)}, arma::vec{v2(i)});
    }
    return result;
  }
  virtual void update_leap_and_shift_indices(arma::uvec& indices, int n_items) = 0;
};

std::unique_ptr<Distance> choose_distance_function(std::string metric);

struct CayleyDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
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
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    indices = arma::regspace<arma::uvec>(0, n_items - 1);
  }
};

struct FootruleDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
    return arma::norm(r1 - r2, 1);
  }
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    return;
  }
};

struct HammingDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
    return arma::sum(r1 != r2);
  }
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    return;
  }
};

struct KendallDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
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
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    return;
  }
};

struct SpearmanDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
    return std::pow(arma::norm(r1 - r2, 2), 2);
  }
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    return;
  }
};

struct UlamDistance : Distance {
  double d(const arma::vec& r1, const arma::vec& r2) override {
    arma::ivec a = arma::conv_to<arma::ivec>::from(r1) - 1;
    arma::ivec b = arma::conv_to<arma::ivec>::from(r2) - 1;
    auto distance = perm0_distance ( a, b );
    return static_cast<double>(distance);
  }
  void update_leap_and_shift_indices(arma::uvec& indices, int n_items) override {
    indices = arma::regspace<arma::uvec>(0, n_items - 1);
  }
};
