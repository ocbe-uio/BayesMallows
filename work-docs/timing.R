library(BayesMallows)
library(tidyverse)
library(latex2exp)

delta_n <- 100
n_items <- ncol(sushi_rankings)
inds <- map(seq_len(nrow(sushi_rankings) / delta_n), function(i) {
  seq(from = (i - 1) * delta_n + 1, to = i * delta_n)
})

smc_mod <- compute_mallows(
  data = setup_rank_data(rankings = sushi_rankings[0, , drop = FALSE]),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000),
  priors = set_priors(lambda = .1)
)

n_particles <- 1000
mcmc_steps <- 5

for(i in 1:10) {
  compute_options <- set_compute_options(
    alpha_prop_sd = .4 - 5e-6 * i * delta_n,
    leap_size = if_else(i < 10, 3, 1))
  priors <- set_priors(lambda = .1)

  smc_mod <- update_mallows(
    smc_mod,
    new_data = setup_rank_data(rankings = sushi_rankings[inds[[i]], ]),
    smc_options = set_smc_options(n_particles = n_particles,
                                  mcmc_steps = mcmc_steps),
    compute_options = compute_options,
    priors = priors
  )
}

rm(smc_mod)
smc_mod <- compute_mallows(
  data = setup_rank_data(rankings = sushi_rankings[0, , drop = FALSE]),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000),
  priors = set_priors(lambda = .1)
)

