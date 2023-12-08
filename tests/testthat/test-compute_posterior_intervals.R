test_that("compute_posterior_intervals works", {
  set.seed(1234)
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 20)
  )
  expect_error(compute_posterior_intervals(mod), "Please specify the burnin.")
  mod$burnin <- 30
  expect_error(compute_posterior_intervals(mod), "burnin < model_fit")

  mod$burnin <- 10
  pi <- compute_posterior_intervals(mod)
  expect_equal(pi$median, "1.669")
  expect_equal(pi$hpdi, "[1.538,1.722]")
  expect_equal(pi$central_interval, "[1.538,1.715]")

  pi <- compute_posterior_intervals(mod, decimals = 1)
  expect_equal(pi$median, "1.7")
  expect_equal(pi$hpdi, "[1.5,1.7]")
  expect_equal(pi$central_interval, "[1.5,1.7]")

  pi <- compute_posterior_intervals(mod, decimals = 2, level = .99)
  expect_equal(pi$median, "1.67")
  expect_equal(pi$hpdi, "[1.54,1.72]")
  expect_equal(pi$central_interval, "[1.54,1.72]")

  expect_error(
    compute_posterior_intervals(mod, parameter = c("alpha", "rho")),
    "'arg' must be of length 1"
  )

  pi <- compute_posterior_intervals(mod, parameter = "rho")
  expect_equal(pi$hpdi[[10]], "[10,11]")
  expect_equal(pi$central_interval[[15]], "[5,6]")

  pi <- compute_posterior_intervals(mod, parameter = "rho", level = .6)
  expect_equal(pi$hpdi[[10]], "[11]")
  expect_equal(pi$central_interval[[15]], "[5,6]")
})

test_that("compute_posterior_intervals works with clusters", {
  set.seed(1234)
  mod <- compute_mallows(
    data = setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 20, burnin = 10)
  )

  pi <- compute_posterior_intervals(mod)
  expect_equal(pi$hpdi[[3]], "[0.831,1.363]")
  expect_equal(pi$central_interval[[2]], "[0.780,1.103]")

  pi <- compute_posterior_intervals(mod, parameter = "rho")
  expect_equal(unique(pi$cluster), factor(paste("Cluster", 1:3)))

  pi <- compute_posterior_intervals(mod,
    parameter = "cluster_probs",
    level = .96, decimals = 4
  )
  expect_equal(pi$hpdi[[2]], "[0.2206,0.4824]")
  expect_equal(pi$central_interval[[1]], "[0.2368,0.4080]")
})

test_that("compute_posterior_intervals works for SMC", {
  set.seed(123)
  mod1 <- compute_mallows(
    setup_rank_data(potato_visual[1:2, ]),
    compute_options = set_compute_options(burnin = 200)
  )
  mod2 <- update_mallows(mod1, new_data = setup_rank_data(potato_visual[3:9, ]))

  pi <- compute_posterior_intervals(mod2)
  expect_equal(pi$parameter, "alpha")
  expect_equal(pi$median, "10.253")
  expect_equal(pi$hpdi, "[8.396,11.858]")
  expect_equal(pi$central_interval, "[8.476,11.985]")

  mod3 <- update_mallows(
    mod2,
    new_data = setup_rank_data(potato_visual[10:12, ])
  )

  pi <- compute_posterior_intervals(mod3, decimals = 2)
  expect_equal(pi$hpdi, "[8.73,11.85]")
  pi <- compute_posterior_intervals(mod3, parameter = "rho")
  expect_equal(pi$hpdi[[20]], "[13,14]")

  expect_error(
    compute_posterior_intervals(mod3, parameter = "cluster_probs"),
    "'arg' should be one of"
  )
})
