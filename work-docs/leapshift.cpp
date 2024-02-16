#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void shift_step(vec& rho_proposal, const vec& rho,
                int u){
  // Shift step:
  double delta_r = rho_proposal(u) - rho(u);
  uvec indices = zeros<uvec>(std::abs(delta_r) + 1);
  indices[0] = u;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      int index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      int index = as_scalar(find(rho == rho(u) + k));
      rho_proposal(index) += 1;
      indices[-(k)] = index;
    }
  }
}

// [[Rcpp::export]]
Rcpp::List leap_and_shift(arma::vec rho, int leap_size){
  vec rho_proposal = rho;
  int n_items = rho.n_elem;

  // Leap 1
  // 1, sample u randomly between 1 and n_items
  ivec a = Rcpp::sample(n_items, 1) - 1;
  int u = a(0);

  // 2, compute the set S for sampling the new rank
  vec support = join_cols(regspace(std::max(1.0, rho(u) - leap_size), 1, rho(u) - 1),
                          regspace(rho(u) + 1, 1, std::min(n_items * 1.0, rho(u) + leap_size)));

  // 3. assign a random element of the support set, this completes the leap step
  ivec b = Rcpp::sample(support.n_elem, 1) - 1;
  int index = b(0);
  // Picked element index-1 from the support set
  rho_proposal(u) = support(index);

  vec support_new = join_cols(regspace(std::max(1.0, rho_proposal(u) - leap_size), 1, rho_proposal(u) - 1),
            regspace(rho_proposal(u) + 1, 1, std::min(n_items * 1.0, rho_proposal(u) + leap_size)));

  // double support_new = std::min(rho_proposal(u) - 1, leap_size * 1.0) +
  //   std::min(n_items - rho_proposal(u), leap_size * 1.0);

  double prob_forward{};
  double prob_backward{};
  if(std::abs(rho_proposal(u) - rho(u)) == 1){
    prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new.size());
    prob_backward = prob_forward;
  } else {
    prob_forward = 1.0 / (n_items * support.n_elem);
    prob_backward = 1.0 / (n_items * support_new.size());
  }
  shift_step(rho_proposal, rho, u);
  return Rcpp::List::create(
    Rcpp::Named("rho_proposal") = rho_proposal,
    Rcpp::Named("prob_forward") = prob_forward,
    Rcpp::Named("prob_backward") = prob_backward,
    Rcpp::Named("u") = u,
    Rcpp::Named("r") = index
  );
}


/*** R
library(tidyverse)
set.seed(1)
n_items <- 6
nmc <- 3000
leap_size <- 3
rho <- seq_len(n_items)
res <- map_dfr(seq_len(nmc), function(i) {
  ls <- leap_and_shift(rho, leap_size)
  tibble(
    distance = sum(abs(ls$rho_proposal - rho)),
    rho_proposal = paste(ls$rho_proposal, collapse = ","),
    prob_forward = ls$prob_forward,
    prob_backward = ls$prob_backward,
    u = ls$u,
    r = ls$r
  )
})

res %>%
  mutate(ur = paste(u, r, sep = ",")) %>%
  group_by(rho_proposal, prob_forward, distance) %>%
  summarise(n_distinct(ur)) %>% View()

res %>%
  group_by(rho_proposal, distance, prob_backward, prob_forward, u, r) %>%
  summarise(pct = n() / nmc) %>%
  View()

res %>%
  filter(rho_proposal == "1,5,4,3,2,6") %>%
  summarise(
    prob_forward = mean(prob_forward),
    prob_backward = mean(prob_backward),
    empirical = n() / nmc
    )

# verify that backward probability is correct
rho_backward <- c(1, 5, 4, 3, 2, 6)
res_backward <- map_dfr(seq_len(nmc), function(i) {
  ls <- leap_and_shift(rho_backward, leap_size)
  tibble(
    rho_proposal = paste(ls$rho_proposal, collapse = ","),
    prob_forward = ls$prob_forward,
    prob_backward = ls$prob_backward
  )
})

res_backward %>%
  filter(rho_proposal == paste(rho, collapse = ",")) %>%
  summarise(
    prob_forward = mean(prob_forward),
    prob_backward = mean(prob_backward),
    empirical = n() / nmc
    )

dd <- res %>%
  group_by(rho_proposal) %>%
  summarise(
    proportion = n() / nmc,
    mean_prob_forward = mean(prob_forward),
    .groups = "drop"
  )

ggplot(dd, aes(x = proportion, y = mean_prob_forward)) +
  geom_point() +
  geom_abline()

ggplot(res, aes(x = prob_forward, y = prob_backward)) +
  geom_point() +
  geom_abline()
*/
