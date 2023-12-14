test_that("expected dist works", {
  expect_equal(
    compute_expected_distance(5, 5, metric = "kendall"), 1.74913673651759
  )
  expect_equal(
    compute_expected_distance(12, 6, metric = "cayley"), 1.3757786382399
  )
  expect_equal(
    compute_expected_distance(1.5 * 7, 7, metric = "hamming"), 2.69246008768975
  )
  expect_equal(
    compute_expected_distance(5 * 30, 30, "ulam"), 4.13353772441732
  )
  expect_equal(
    compute_expected_distance(3.5 * 45, 45, "footrule"), 0.0804585327652995
  )
  expect_equal(
    compute_expected_distance(4 * 10, 10, "spearman"), 0.00603271004508636
  )
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
    "n_items must be a positive integer"
  )
  expect_error(
    compute_expected_distance(-2, 15, "ulam"),
    "alpha must be a non-negative value"
  )
})
