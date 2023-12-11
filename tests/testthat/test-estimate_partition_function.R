test_that("estimate_partition_function works", {
  set.seed(1)
  alpha_vector <- seq(from = 0, to = 10, by = 0.5)
  n_items <- 20

  metrics <- c("footrule", "spearman", "kendall", "cayley", "hamming", "ulam")
  expectations <- c(-0.196249967926383, 16.7714664134817,
                    -13.9827539118596, 15.4312563640867,
                    -30.5913961068133, -3.21589663450721)
  names(expectations) <- metrics

  for(m in metrics) {
    fit <- estimate_partition_function(
      method = "importance_sampling",
      alpha_vector = alpha_vector,
      n_items = n_items,
      metric = m,
      n_iterations = 1e3)
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
