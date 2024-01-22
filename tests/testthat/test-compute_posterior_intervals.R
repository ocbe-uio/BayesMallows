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
  expect_equal(pi$median, "0.861")
  expect_equal(pi$hpdi, "[0.814,1.083]")
  expect_equal(pi$central_interval, "[0.819,1.080]")

  pi <- compute_posterior_intervals(mod, decimals = 1)
  expect_equal(pi$median, "0.9")
  expect_equal(pi$hpdi, "[0.8,1.1]")
  expect_equal(pi$central_interval, "[0.8,1.1]")

  pi <- compute_posterior_intervals(mod, decimals = 2, level = .99)
  expect_equal(pi$median, "0.86")
  expect_equal(pi$hpdi, "[0.81,1.08]")
  expect_equal(pi$central_interval, "[0.81,1.08]")

  expect_error(
    compute_posterior_intervals(mod, parameter = c("alpha", "rho")),
    "'arg' must be of length 1"
  )

  pi <- compute_posterior_intervals(mod, parameter = "rho")
  expect_equal(pi$hpdi[[10]], "[6,7]")
  expect_equal(pi$central_interval[[15]], "[1,2]")

  pi <- compute_posterior_intervals(mod, parameter = "rho", level = .6)
  expect_equal(pi$hpdi[[10]], "[7]")
  expect_equal(pi$central_interval[[15]], "[1,2]")
})

test_that("compute_posterior_intervals works with clusters", {
  set.seed(1234)
  mod <- compute_mallows(
    data = setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 20, burnin = 10)
  )

  pi <- compute_posterior_intervals(mod)
  expect_equal(pi$hpdi[[3]], "[0.608,0.753]")
  expect_equal(pi$central_interval[[2]], "[0.710,1.140]")

  pi <- compute_posterior_intervals(mod, parameter = "rho")
  expect_equal(unique(pi$cluster), factor(paste("Cluster", 1:3)))

  pi <- compute_posterior_intervals(mod,
    parameter = "cluster_probs",
    level = .96, decimals = 4
  )
  expect_equal(pi$hpdi[[2]], "[0.2552,0.3708]")
  expect_equal(pi$central_interval[[1]], "[0.3184,0.4840]")
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
  expect_equal(pi$median, "2.744")
  expect_equal(pi$hpdi, "[2.015,3.476]")
  expect_equal(pi$central_interval, "[2.047,3.530]")

  mod3 <- update_mallows(
    mod2,
    new_data = setup_rank_data(potato_visual[10:12, ])
  )

  pi <- compute_posterior_intervals(mod3, decimals = 2)
  expect_equal(pi$hpdi, "[2.18,3.75]")
  pi <- compute_posterior_intervals(mod3, parameter = "rho")
  expect_equal(pi$hpdi[[20]], "[1,5]")

  expect_error(
    compute_posterior_intervals(mod3, parameter = "cluster_probs"),
    "'arg' should be one of"
  )
})
