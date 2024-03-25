library(BayesMallows)
m0 <- compute_mallows(
  data = setup_rank_data(preferences = subset(beach_preferences, assessor == 2)),
  compute_options = set_compute_options(nmc = 2000, burnin = 500)
)

m1 <- update_mallows(
  model = m0,
  new_data = setup_rank_data(
    preferences = beach_preferences[1:10, ],
    user_ids = 1,
    n_items = 15
  ),
  smc_options = set_smc_options(n_particles = 2)
)

m2 <- update_mallows(
  model = m1,
  new_data = setup_rank_data(
    preferences = beach_preferences[11:20, ],
    user_ids = 1
  ),
  smc_options = set_smc_options(n_particles = 2)
)
