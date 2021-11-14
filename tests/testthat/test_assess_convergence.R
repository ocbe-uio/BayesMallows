context("Testing assess_convergence")

test_that("assess_convergence fails when it should", {
  m <- compute_mallows(potato_visual, nmc = 10)
  class(m) <- NULL
  expect_error(assess_convergence(m))
  expect_error(assess_convergence(m, parameter = "aaa"))
  expect_error(assess_convergence(m, parameter = "Rtilde"))
  m <- compute_mallows(potato_visual, nmc = 10, save_aug = TRUE)
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = 1:100))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = c("a", "b")))
  expect_s3_class(assess_convergence(m), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "alpha"), "ggplot")

  expect_s3_class(assess_convergence(m, parameter = "rho"), "ggplot")

  m <- compute_mallows(preferences = beach_preferences, nmc = 10, error_model = "bernoulli")
  expect_s3_class(assess_convergence(m, parameter = "theta"), "ggplot")
})


test_that("assess_convergence works with mixtures", {
  m <- compute_mallows_mixtures(n_clusters = 3, rankings = potato_visual, nmc = 10)
  class(m) <- NULL
  expect_error(assess_convergence(m))
  expect_error(assess_convergence(m, parameter = "aaa"))
  expect_error(assess_convergence(m, parameter = "Rtilde"))

  m <- compute_mallows_mixtures(n_clusters = 3, rankings = potato_visual, nmc = 10,
                                save_aug = TRUE)
  expect_s3_class(assess_convergence(m, parameter = "alpha"), "ggplot")
  expect_error(assess_convergence(m, parameter = "Rtilde"))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = 1:100))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = c("a", "b")))

  expect_s3_class(assess_convergence(m), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "rho"), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "cluster_probs"), "ggplot")

  m <- compute_mallows(potato_visual, nmc = 10, n_clusters = 2, save_clus = TRUE)
  expect_s3_class(assess_convergence(m, parameter = "cluster_probs"), "ggplot")
})
