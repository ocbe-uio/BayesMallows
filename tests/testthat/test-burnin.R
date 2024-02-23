test_that("burnin works", {
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 30, burnin = 20)
  )
  expect_equal(burnin(mod), 20)
  burnin(mod) <- 10
  expect_equal(burnin(mod), 10)
  expect_error(burnin(mod) <- 40, "Burnin cannot be larger")

  mod <- update_mallows(mod, new_data = setup_rank_data(potato_weighing))
  expect_equal(burnin(mod), 0)
  expect_error(burnin(mod) <- 3, "Cannot set burnin for SMC model.")

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 30)
  )
  expect_equal(burnin(mod), NULL)

  mod <- compute_mallows_mixtures(
    n_clusters = 1:3,
    data = setup_rank_data(cluster_data),
    compute_options = set_compute_options(nmc = 20)
  )
  expect_equal(burnin(mod), list(NULL, NULL, NULL))
  burnin(mod) <- 3
  expect_equal(burnin(mod), list(3, 3, 3))
  burnin(mod) <- 4:6
  expect_equal(burnin(mod), list(4, 5, 6))

})
