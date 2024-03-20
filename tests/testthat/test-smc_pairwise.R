test_that("update_mallows works with pairwise preferences", {
  set.seed(3)
  dat <- subset(beach_preferences, assessor <= 10)
  mod_init <- compute_mallows(
    data = setup_rank_data(
      preferences = beach_preferences
    ),
    compute_options = set_compute_options(nmc = 3000, burnin = 1000)
  )

  # Next we provide assessors 21 to 60 one at a time.
  mod <- mod_init
  for (i in 21:22) {
    mod <- update_mallows(
      model = mod,
      new_data = setup_rank_data(
        preferences = subset(beach_preferences, assessor == i),
        user_ids = i
      ),
      smc_options = set_smc_options(latent_sampling_lag = 0)
    )
  }

  expect_equal(mean(mod$alpha_samples), 4.12374185445437)
  expect_equal(sd(mod$alpha_samples), 0.858531920995246)

  mod <- mod_init
  for (i in 23:24) {
    mod <- update_mallows(
      model = mod,
      new_data = setup_rank_data(
        preferences = subset(beach_preferences, assessor == i),
        user_ids = i
      ),
      smc_options = set_smc_options(latent_sampling_lag = 1, max_topological_sorts = 10)
    )
  }

  expect_equal(mean(mod$alpha_samples), 2.89026770205348)
  expect_equal(sd(mod$alpha_samples), 0.648145423252438)
})
