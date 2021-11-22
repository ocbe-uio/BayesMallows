context("Testing print.BayesMallows and print.BayesMallowsMixtures")

test_that("print.BayesMallows fails when it should", {
  class(mtcars) <- "BayesMallows"
  expect_error(print(mtcars))
  m <- compute_mallows(potato_visual, nmc = 5)
  m$n_items <- NULL
  expect_error(print(m), "BayesMallows object must have elements n_items and n_assessors.")
  m$n_assessors <- NULL
  expect_error(print(m),
               "BayesMallows object must have elements n_items and n_assessors.")
})


test_that("print.BayesMallowsMixtures fails when it should", {
  class(mtcars) <- "BayesMallowsMixtures"
  expect_error(print(mtcars),
               regexp = "All elements of a BayesMallowsMixtures object must be of class BayesMallows.")
  m <- compute_mallows_mixtures(potato_visual, n_clusters = 1:3, nmc = 5)
  class(m[[1]]) <- "list"
  expect_error(print(m),
               regexp = "All elements of a BayesMallowsMixtures object must be of class BayesMallows.")

  m <- compute_mallows_mixtures(potato_visual, n_clusters = 1:3, nmc = 5)
  expect_output(print(m),
                "Collection of 3 Bayesian Mallows Mixture Models with the following number of mixture components:")
})
