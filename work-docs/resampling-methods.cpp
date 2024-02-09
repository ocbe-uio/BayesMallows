#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::ivec count_to_index(const arma::vec& counts) {
  arma::ivec indices(arma::sum(counts));

  int idx = 0;
  for (arma::uword i = 0; i < counts.n_elem; ++i) {
    for (int j = 0; j < counts(i); ++j) {
      indices(idx++) = i;
    }
  }
  return indices;
}

arma::ivec digitize(const arma::vec& bins, const arma::vec& values) {
  arma::ivec indices(values.n_elem);

  for (arma::uword i = 0; i < values.n_elem; ++i) {
    double val = values(i);
    arma::uword index = 0;
    while (index < bins.n_elem && val >= bins(index)) {
      ++index;
    }
    indices(i) = index;
  }

  return indices;
}

arma::ivec stratsys(arma::vec probs, bool stratified) {
  size_t N = probs.size();
  arma::vec u(N);
  Rcpp::NumericVector rn;
  if(stratified) {
    rn = Rcpp::runif(N, 0, 1);
  } else {
    rn = Rcpp::NumericVector(N, R::runif(0, 1));
  }

  for(size_t i{}; i < N; i++) u(i) = (i + rn(i)) / N;
  return digitize(arma::cumsum(probs), u);
}

// [[Rcpp::export]]
arma::ivec multinomial_resampling(arma::vec probs) {
  return Rcpp::sample(
    probs.size(), probs.size(), true,
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(probs)), false);
}

// [[Rcpp::export]]
arma::ivec residual_resampling(arma::vec probs) {
  int N = probs.size();
  arma::vec counts = arma::floor(probs * N);
  int R = N - arma::sum(counts);
  arma::ivec result = count_to_index(counts);

  arma::vec residual_weights = probs - counts / N;
  residual_weights = residual_weights / arma::sum(residual_weights);
  arma::ivec part2 = Rcpp::sample(
    N, R, true,
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(residual_weights)), false);

  return arma::join_cols(result, part2);
}

// [[Rcpp::export]]
arma::ivec stratified_resampling(arma::vec probs) {
  return stratsys(probs, true);
}

// [[Rcpp::export]]
arma::ivec systematic_resampling(arma::vec probs) {
  return stratsys(probs, false);
}


/*** R
library(tidyverse)

n <- 4000
probs <- rexp(n, rate = .0001)
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
