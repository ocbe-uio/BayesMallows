test_that("compute_exact_partition_function works", {
  expect_error(
    compute_exact_partition_function(3, -1),
    "n_items must be a positive integer"
  )
  expect_error(
    compute_exact_partition_function(3, 2.3),
    "n_items must be a positive integer"
  )
  expect_error(
    compute_exact_partition_function(3, 0),
    "n_items must be a strictly positive number of length one"
  )
  expect_error(
    compute_exact_partition_function(-2, 3),
    "alpha must be a strictly positive number of length one"
  )
  expect_error(
    compute_exact_partition_function(rnorm(2), 3),
    "alpha must be a strictly positive number of length one"
  )

  expect_equal(
    compute_exact_partition_function(alpha = 3, n_items = 10, metric = "cayley"),
    13.0481794289176
  )
  expect_equal(
    compute_exact_partition_function(alpha = 3, n_items = 10, metric = "hamming"),
    12.4542713806513
  )
  expect_equal(
    compute_exact_partition_function(alpha = 3, n_items = 10, metric = "kendall"),
    9.69641008390133
  )
})
