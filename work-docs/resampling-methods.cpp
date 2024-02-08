#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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

// [[Rcpp::export]]
Rcpp::IntegerVector multinomial_resampling(arma::vec probs) {
  return Rcpp::sample(probs.size(), probs.size(), true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(probs)), false);
}

// [[Rcpp::export]]
Rcpp::IntegerVector residual_resampling(arma::vec probs) {
  int N = probs.size();
  // Stage 1, deterministic
  arma::vec counts = arma::floor(probs * N);
  int R = N - arma::sum(counts);
  Rcpp::IntegerVector result = count_to_index(counts);

  // Stage 2, random
  arma::vec residual_weights = probs - counts / N;
  residual_weights = residual_weights / arma::sum(residual_weights);
  Rcpp::IntegerVector part2 = Rcpp::sample(N, R, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(residual_weights)), false);
  for(auto x : part2) result.push_back(x);

  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector stratified_resampling(arma::vec probs, bool systematic = false) {
  size_t N = probs.size();
  arma::vec u(N);
  if(systematic) {
    double rn = R::runif(0, 1);
    for(size_t i{}; i < N; i++) {
      u(i) = (i + rn) / N;
    }
  } else {
    for(size_t i{}; i < N; i++) {
      u(i) = (i + R::runif(0, 1)) / N;
    }
  }

  arma::vec bins = arma::join_cols(arma::vec{0}, arma::cumsum(probs));
  arma::uvec counts = arma::histc(u, bins);
  Rcpp::IntegerVector result = count_to_index(counts);

  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector systematic_resampling(arma::vec probs) {
  return stratified_resampling(probs, true);
}


/*** R
library(tidyverse)

n <- 3000
probs <- rexp(n, rate = .1)
probs <- probs / sum(probs)

samplefun <- function(f) {
  as.numeric(table(factor(f(probs), levels = 0:(n-1))))
}

dat <- tibble(
  expected_value = n * probs,
  multinomial_resampling = samplefun(multinomial_resampling),
  residual_resampling = samplefun(residual_resampling),
  stratified_resampling = samplefun(stratified_resampling),
  systematic_resampling = samplefun(systematic_resampling)
)

dat %>%
  pivot_longer(cols = -expected_value) %>%
  ggplot(aes(x = value, y = expected_value)) +
  geom_point() +
  geom_abline() +
  facet_wrap(vars(name))
*/
