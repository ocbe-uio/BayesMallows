test_that("print function works", {
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual[1:5, ]),
    compute_options = set_compute_options(nmc = 100, burnin = 0)
  )

  expect_snapshot(mod)

  mod2 <- update_mallows(
    model = mod,
    new_data = setup_rank_data(potato_visual[6:12, ]),
    smc_options = set_smc_options(n_particles = 3)
  )

  expect_snapshot(mod2)

  mod3 <- compute_mallows_mixtures(
    n_clusters = 1:3,
    data = setup_rank_data(cluster_data),
    compute_options = set_compute_options(nmc = 30, burnin = 0)
  )

  expect_snapshot(mod3)
})
