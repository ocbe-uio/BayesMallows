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

  expect_equal(mean(mod$alpha_samples), 3.99645235790935)
  expect_equal(sd(mod$alpha_samples), 0.872720087078698)

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

  expect_equal(mean(mod$alpha_samples), 3.03458072249661)
  expect_equal(sd(mod$alpha_samples), 0.670792690537946)
})

test_that("update_mallows works with existing users updating their data", {
  set.seed(22)
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
    )
  )

  expect_equal(
    mean(m1$alpha_samples),
    2.75891148770247
  )

  expect_equal(
    order(apply(m1$augmented_rankings, 1, mean)),
    c(1L, 6L, 9L, 3L, 7L, 4L, 10L, 11L, 8L, 12L, 15L, 13L, 14L, 2L, 5L)
  )

  m2 <- update_mallows(
    model = m1,
    new_data = setup_rank_data(
      preferences = beach_preferences[11:20, ],
      user_ids = 1
    )
  )

  expect_equal(
    mean(m2$alpha_samples),
    2.82143649083478
  )

  expect_equal(
    order(apply(m2$augmented_rankings, 1, mean)),
    c(6L, 4L, 10L, 1L, 2L, 7L, 8L, 11L, 5L, 3L, 12L, 9L, 13L, 15L, 14L)
  )
})

test_that("update_mallows works with both new and updated users", {
  set.seed(22)
  m0 <- compute_mallows(
    data = setup_rank_data(preferences = subset(beach_preferences, assessor == 2)),
    compute_options = set_compute_options(nmc = 2000, burnin = 500)
  )

  m1 <- update_mallows(
    model = m0,
    new_data = setup_rank_data(
      preferences = beach_preferences[c(1:10, 51:60), ],
      user_ids = c(1, 3),
      n_items = 15
    )
  )

  m2 <- update_mallows(
    model = m1,
    new_data = setup_rank_data(
      preferences = subset(beach_preferences, assessor %in% c(1, 3, 4)),
      user_ids = c(1, 3, 4)
    )
  )

  expect_equal(m1$data$user_ids, c(1, 3))
  expect_equal(m2$data$user_ids, c(1, 3, 4))
})
