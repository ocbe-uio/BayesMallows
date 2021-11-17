test_that("print function works", {
  m <- compute_mallows(potato_visual, nmc = 3)
  expect_output(print(m), "Bayesian Mallows Model with 20 items and 12 assessors.")
})
