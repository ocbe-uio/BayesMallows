context("Testing plot.BayesMallows")

test_that("plot.BayesMallows fails when it should", {
  class(mtcars) <- "BayesMallows"
  expect_error(plot(mtcars))
  m <- compute_mallows(potato_visual, nmc = 100)
  expect_error(plot(m))
  expect_error(plot(m, parameter = "Rtilde", burnin = 50))
  expect_error(plot(m, burnin = 1e7))
})


test_that("plot.BayesMallows works", {
  m <- compute_mallows(potato_visual, nmc = 10)
  expect_s3_class(plot(m, burnin = 3), "ggplot")
  m$burnin <- 4
  expect_s3_class(plot(m), "ggplot")
  expect_s3_class(plot(m, parameter = "rho"), "ggplot")
  expect_s3_class(plot(m, parameter = "rho", items = 2:3), "ggplot")

  m <- compute_mallows(potato_visual, nmc = 10, n_clusters = 3, save_clus = TRUE)
  expect_s3_class(plot(m, burnin = 4, parameter = "cluster_probs"), "ggplot")
  expect_s3_class(plot(m, burnin = 4, parameter = "cluster_assignment"), "ggplot")

  m <- compute_mallows(preferences = beach_preferences[1:100, ], nmc = 10, error_model = "bernoulli")
  expect_s3_class(plot(m, burnin = 3, parameter = "theta"), "ggplot")
})
