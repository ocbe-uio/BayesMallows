test_that("estimate_partition_function works", {
  set.seed(1)
  alpha_vector <- seq(from = 0, to = 10, by = 0.5)
  n_items <- 20

  metrics <- c("footrule", "spearman", "kendall", "cayley", "hamming", "ulam")
  expectations <- c(
    0.4961490378154, 19.75045734511,
    22.3471310441124, -9.13330233602348,
    -78.559587447592, -10.9508829949516
  )
  names(expectations) <- metrics

  for (m in metrics) {
    fit <- estimate_partition_function(
      method = "importance_sampling",
      alpha_vector = alpha_vector,
      n_items = n_items,
      metric = m,
      n_iterations = 1e3
    )
    expect_equal(fit[[5]], expectations[[m]])
  }

  fit <- estimate_partition_function(
    method = "asymptotic",
    alpha_vector = alpha_vector,
    n_items = n_items,
    metric = "footrule",
    n_iterations = 50
  )
  expect_equal(fit[[4]], 0.0151144482198799)

  fit <- estimate_partition_function(
    method = "asymptotic",
    alpha_vector = alpha_vector,
    n_items = n_items,
    metric = "spearman",
    n_iterations = 50
  )
  expect_equal(fit[[4]], -41.9129447085325)
})

test_that("estimate_partition_function works in parallel", {
  cl <- parallel::makeCluster(2)
  set.seed(1)
  alpha_vector <- seq(from = 0, to = 10, by = 0.5)

  fit <- estimate_partition_function(
    method = "importance_sampling",
    alpha_vector = alpha_vector,
    n_items = 34,
    metric = "kendall",
    n_iterations = 1e3,
    cl = cl
  )
  expect_equal(fit[[3]], 33.4259948749862)
  parallel::stopCluster(cl)

  mod <- compute_mallows(
    data = setup_rank_data(t(replicate(3, sample(34)))),
    model_options = set_model_options(metric = "kendall"),
    compute_options = set_compute_options(nmc = 10),
    logz_estimate = fit
  )
  expect_s3_class(mod, "BayesMallows")
})
