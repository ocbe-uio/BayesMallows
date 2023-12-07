cl <- parallel::makeCluster(2)
test_that("compute_mallows works with seed in parallel", {
  set.seed(1)
  mod <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10),
    cl = cl
  )
  expect_equal(mod$rho$value[[145]], 4)
  expect_equal(mod$rho$value[[89]], 3)
})

test_that("compute_mallows works with initial values and clusters", {
  set.seed(1)
  mod <- compute_mallows(
    setup_rank_data(potato_visual),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 10),
    initial_values = set_initial_values(alpha_init = c(2, 3)),
    cl = cl
  )
  expect_true(
    all(subset(mod$alpha, chain == 1 & iteration == 1)$value == 2)
  )
  expect_true(
    all(subset(mod$alpha, chain == 2 & iteration == 1)$value == 3)
  )
})
parallel::stopCluster(cl)
