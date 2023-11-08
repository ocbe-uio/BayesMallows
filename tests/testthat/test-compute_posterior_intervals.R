test_that("compute posterior intervals works", {
  set.seed(123)
  m <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10)
  )
  expect_error(compute_posterior_intervals(m))
  expect_error(compute_posterior_intervals(m, burnin = 100))
  expect_error(compute_posterior_intervals(m, burnin = 7, parameter = "dsdsd"))
  expect_error(compute_posterior_intervals(m, burnin = 7, level = 1.2))
  expect_snapshot(
    compute_posterior_intervals(m, burnin = 7, level = .05, parameter = "alpha")
  )

  expect_snapshot(
    compute_posterior_intervals(m, burnin = 7, level = .1, parameter = "cluster_probs")
  )

  expect_snapshot(
    compute_posterior_intervals(m, burnin = 7, level = .01, parameter = "rho")
  )


  set.seed(22)
  m <- compute_mallows(
    setup_rank_data(potato_visual),
    model = set_model_options(n_clusters = 2),
    compute_options = set_compute_options(nmc = 10)
  )
  expect_snapshot(
    compute_posterior_intervals(m, burnin = 8)
  )
})
