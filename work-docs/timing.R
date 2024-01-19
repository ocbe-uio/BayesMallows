library(BayesMallows)
set.seed(123)

mod_bmm <- compute_mallows(
  data = setup_rank_data(sushi_rankings),
  compute_options = set_compute_options(nmc = 2000, burnin = 200)
)

mod_init <- compute_mallows(
  data = setup_rank_data(sushi_rankings[1:100, ]),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

system.time({
  mod_smc <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = sushi_rankings[101:2000, ]),
    smc_options = set_smc_options(n_particles = 2000, mcmc_steps = 5)
  )
})
# 0.25 in parallel
# 0.31 sequentially
system.time({
  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(sushi_rankings[2001:5000, ])
  )
})
# 0.77 in parallel
# 0.82 sequentially
