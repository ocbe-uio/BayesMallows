test_that("acceptance rates are computed", {
  library(parallel)
  cl <- makeCluster(2)
  set.seed(1)
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100),
    cl = cl
  )
  expect_equal(mod$alpha_acceptance, list(0.911111111111111, 0.833333333333333))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 100),
    cl = cl
  )
  expect_equal(mod$alpha_acceptance, list(0.0444444444444444, 0.0333333333333333))

  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(burnin = 10, nmc = 100, alpha_prop_sd = 0.00001),
    cl = cl
  )
  expect_equal(mod$alpha_acceptance, list(0.988888888888889, 0.988888888888889))
  stopCluster(cl)
})
