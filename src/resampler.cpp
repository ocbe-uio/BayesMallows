#include "resampler.h"

template<typename T>
Rcpp::IntegerVector count_to_index(const T& counts) {
  Rcpp::IntegerVector result;
  for(int i{}; i < counts.size(); i++) {
    for(int j{}; j < counts[i]; j++) {
      result.push_back(i);
    }
  }
  return result;
}

Rcpp::IntegerVector Multinomial::resample(arma::vec probs) {
  return Rcpp::sample(
    probs.size(), probs.size(), true,
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(probs)), false);
}

Rcpp::IntegerVector Residual::resample(arma::vec probs) {
  int N = probs.size();
  arma::vec counts = arma::floor(probs * N);
  int R = N - arma::sum(counts);
  Rcpp::IntegerVector result = count_to_index(counts);

  arma::vec residual_weights = probs - counts / N;
  residual_weights = residual_weights / arma::sum(residual_weights);
  Rcpp::IntegerVector part2 = Rcpp::sample(N, R, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(residual_weights)), false);
  for(auto x : part2) result.push_back(x);

  return result;
}

Rcpp::IntegerVector stratsys(arma::vec probs, bool stratified) {
  size_t N = probs.size();
  arma::vec u(N);
  Rcpp::NumericVector rn;
  if(stratified) {
    rn = Rcpp::runif(N, 0, 1);
  } else {
    rn = Rcpp::NumericVector(N, R::runif(0, 1));
  }
  for(size_t i{}; i < N; i++) u(i) = (i + rn(i)) / N;

  arma::vec bins = arma::join_cols(arma::vec{0}, arma::cumsum(probs));
  arma::uvec counts = arma::histc(u, bins);
  Rcpp::IntegerVector result = count_to_index(counts);

  return result;
}


Rcpp::IntegerVector Stratified::resample(arma::vec probs) {
  return stratsys(probs, true);
}

Rcpp::IntegerVector Systematic::resample(arma::vec probs) {
  return stratsys(probs, false);
}
