context("Testing rank_conversion")

test_that("create_ranking returns correct values", {
  orderings <- matrix(c(1, 2, 2, 3, 3, 1), ncol = 3)
  rankings <- matrix(c(1, 3, 2, 1, 3, 2), ncol = 3)
  expect_equal(create_ranking(orderings), rankings)
  orderings[2, 2] <- NA_real_
  rankings[2, 3] <- NA_real_
  expect_equal(create_ranking(orderings), rankings)
  orderings[2, 3] <- NA_real_
  rankings[2, 1] <- NA_real_
  expect_equal(create_ranking(orderings), rankings)
})


test_that("create_orderings returns correct values", {
  orderings <- matrix(c(1, 2, 2, 3, 3, 1), ncol = 3)
  rankings <- matrix(c(1, 3, 2, 1, 3, 2), ncol = 3)
  expect_equal(create_ordering(rankings), orderings)
  orderings[2, 2] <- NA_real_
  rankings[2, 3] <- NA_real_
  expect_equal(create_ordering(rankings), orderings)
  orderings[2, 3] <- NA_real_
  rankings[2, 1] <- NA_real_
  expect_equal(create_ordering(rankings), orderings)
})
