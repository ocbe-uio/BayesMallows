test_that("acceptance rates are computed", {
  cl <- parallel::makeCluster(2)
  set.seed(1)
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100),
    cl = cl
  )
  expect_equal(get_acceptance_ratios(mod)$alpha_acceptance, list(0.911111111111111, 0.833333333333333))
  expect_equal(get_acceptance_ratios(mod)$rho_acceptance, list(0.977777777777778, 0.9))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 100),
    cl = cl
  )
  expect_equal(get_acceptance_ratios(mod)$alpha_acceptance, list(0.0444444444444444, 0.0333333333333333))
  expect_equal(get_acceptance_ratios(mod)$rho_acceptance, list(0.977777777777778, 0.988888888888889))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 0.00001),
    cl = cl
  )
  expect_equal(get_acceptance_ratios(mod)$alpha_acceptance, list(0.988888888888889, 0.988888888888889))
  expect_equal(get_acceptance_ratios(mod)$rho_acceptance, list(0.922222222222222, 0.888888888888889))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, leap_size = 9),
    cl = cl
  )
  expect_equal(get_acceptance_ratios(mod)$alpha_acceptance, list(0.711111111111111, 0.777777777777778))
  expect_equal(get_acceptance_ratios(mod)$rho_acceptance, list(0.566666666666667, 0.511111111111111))

  set.seed(1)
  mod1 <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, rho_thinning = 1),
    cl = cl
  )
  set.seed(1)
  mod2 <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, rho_thinning = 3),
    cl = cl
  )
  expect_equal(mod1$acceptance_ratios$rho_acceptance, mod2$acceptance_ratios$rho_acceptance)
  parallel::stopCluster(cl)
})

test_that("acceptance rates work with SMC", {
  set.seed(1)
  mod0 <- compute_mallows(
    data = setup_rank_data(potato_visual[1:6, ]),
    compute_options = set_compute_options(nmc = 100, burnin = 10)
  )

  mod1 <- update_mallows(mod0, setup_rank_data(potato_visual[7:12, ]),
    smc_options = set_smc_options(n_particles = 100)
  )
  expect_equal(as.numeric(mod1$acceptance_ratios$alpha_acceptance), .966)
  expect_equal(as.numeric(mod1$acceptance_ratios$rho_acceptance), .964)

  mod2 <- update_mallows(mod1, setup_rank_data(potato_weighing[1:6, ]))
  expect_equal(as.numeric(mod2$acceptance_ratios$alpha_acceptance), .958)
  expect_equal(as.numeric(mod2$acceptance_ratios$rho_acceptance), .888)

  mod3 <- update_mallows(mod2, setup_rank_data(potato_weighing[7:12, ]))
  expect_equal(as.numeric(mod3$acceptance_ratios$alpha_acceptance), .902)
  expect_equal(as.numeric(mod3$acceptance_ratios$rho_acceptance), .854)
})
