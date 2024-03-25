library(BayesMallows)
m0 <- compute_mallows(
  data = setup_rank_data(preferences = subset(beach_preferences, assessor == 2 & bottom_item < 5 & top_item < 5),
                         n_items = 5),
  compute_options = set_compute_options(nmc = 2000, burnin = 500)
)

m1 <- update_mallows(
  model = m0,
  new_data = setup_rank_data(
    preferences = data.frame(
      assessor = 1,
      bottom_item = c(5, 4, 3, 2),
      top_item = c(4, 3, 2, 1)
    ),
    user_ids = 1,
    n_items = 5
  ),
  smc_options = set_smc_options(n_particles = 2)
)

m2 <- update_mallows(
  model = m1,
  setup_rank_data(
    preferences = data.frame(
      assessor = 1,
      bottom_item = c(5, 4, 3, 2, 5),
      top_item = c(4, 3, 2, 1, 3)
    ),
    user_ids = 1,
    n_items = 5
  ),
  smc_options = set_smc_options(n_particles = 2)
)
