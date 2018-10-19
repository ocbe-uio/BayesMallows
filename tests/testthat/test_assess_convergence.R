context("Testing assess_convergence")

test_that("assess_convergence fails when it should", {
  m <- compute_mallows(potato_visual, nmc = 10)
  class(m) <- NULL
  expect_error(assess_convergence(m))
  m <- compute_mallows(potato_visual, nmc = 10)
  expect_error(assess_convergence(m, type = "aaa"))
  m <- compute_mallows(potato_visual, nmc = 10)
  expect_error(assess_convergence(m, type = "Rtilde"))
  m <- compute_mallows(potato_visual, nmc = 10, save_aug = TRUE)
  expect_error(assess_convergence(m, type = "Rtilde", assessors = 1:100))
  expect_error(assess_convergence(m, type = "Rtilde", assessors = c("a", "b")))
})
