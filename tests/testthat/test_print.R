context("Testing print.BayesMallows and print.BayesMallowsMixtures")

test_that("print.BayesMallows fails when it should", {
  class(mtcars) <- "BayesMallows"
  expect_error(print(mtcars))
  m <- compute_mallows(potato_visual, nmc = 100)
  m$n_items <- NULL
  expect_error(print(m))
  m$n_assessors <- NULL
  expect_error(print(m))
})


test_that("print.BayesMallowsMixtures fails when it should", {
  class(mtcars) <- "BayesMallowsMixtures"
  expect_error(print(mtcars))
  m <- compute_mallows_mixtures(potato_visual, n_clusters = 1:3, nmc = 100)
  class(m[[1]]) <- "list"
  expect_error(print(m))
})
