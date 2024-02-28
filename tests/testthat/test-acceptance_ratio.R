test_that("acceptance rates are computed", {
  library(parallel)
  cl <- makeCluster(2)
  set.seed(1)
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100),
    cl = cl
  )
  expect_equal(mod$acceptance_ratios$alpha_acceptance, list(0.911111111111111, 0.833333333333333))
  expect_equal(mod$acceptance_ratios$rho_acceptance, list(0.977777777777778, 0.9))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 100),
    cl = cl
  )
  expect_equal(mod$acceptance_ratios$alpha_acceptance, list(0.0444444444444444, 0.0333333333333333))
  expect_equal(mod$acceptance_ratios$rho_acceptance, list(0.977777777777778, 0.988888888888889))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 0.00001),
    cl = cl
  )
  expect_equal(mod$acceptance_ratios$alpha_acceptance, list(0.988888888888889, 0.988888888888889))
  expect_equal(mod$acceptance_ratios$rho_acceptance, list(0.922222222222222, 0.888888888888889))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, leap_size = 9),
    cl = cl
  )
  expect_equal(mod$acceptance_ratios$alpha_acceptance, list(0.711111111111111, 0.777777777777778))
  expect_equal(mod$acceptance_ratios$rho_acceptance, list(0.566666666666667, 0.511111111111111))

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
  stopCluster(cl)
})
