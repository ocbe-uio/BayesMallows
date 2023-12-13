test_that("expected dist works", {
  expect_equal(round(compute_expected_distance(5, 5, metric = "kendall"), 6), 1.749137)
  expect_equal(round(compute_expected_distance(12, 6, metric = "cayley"), 6), 1.375779)
  expect_equal(round(compute_expected_distance(1.5 * 7, 7, metric = "hamming"), 6), 2.69246)
  expect_equal(round(compute_expected_distance(5 * 30, 30, "ulam"), 6), 4.133538)
  expect_equal(round(compute_expected_distance(3.5 * 45, 45, "footrule"), 6), 0.080459)
  expect_equal(round(compute_expected_distance(4 * 10, 10, "spearman"), 6), 0.006033)
})


test_that("expected dist fails when it should", {
  expect_error(
    compute_expected_distance(10, 150, "footrule"),
    "Not available for requested number of items."
  )
  expect_error(
    compute_expected_distance(10, 150, "ulam"),
    "Not available for requested number of items."
  )
  expect_error(
    compute_expected_distance(10, -2, "spearman"),
    "Number of items must be a positive integer"
  )
  expect_error(
    compute_expected_distance(-2, 15, "ulam"),
    "alpha must be a non-negative value"
  )
})
