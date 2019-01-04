context("Testing plot.BayesMallows")

test_that("plot.BayesMallows fails when it should", {
  class(mtcars) <- "BayesMallows"
  expect_error(plot(mtcars))
  m <- compute_mallows(potato_visual, nmc = 100)
  expect_error(plot(m))
  expect_error(plot(m, parameter = "Rtilde", burnin = 50))
  expect_error(plot(m, burnin = 1e7))
})
