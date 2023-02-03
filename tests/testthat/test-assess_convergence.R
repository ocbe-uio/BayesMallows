context("Testing assess_convergence")

test_that("assess_convergence fails when it should", {
  m <- compute_mallows(potato_visual, nmc = 10)
  class(m) <- NULL
  expect_error(assess_convergence(m))
  expect_error(assess_convergence(m, parameter = "aaa"))
  expect_error(assess_convergence(m, parameter = "Rtilde"))
  m <- compute_mallows(potato_visual, nmc = 10, save_aug = TRUE)
  expect_error(assess_convergence(m, parameter = "theta"), regexp = "Theta not available")
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = 1:100))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = c("a", "b")))
  expect_s3_class(assess_convergence(m), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "alpha"), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "rho"), "ggplot")
  expect_message(
    assess_convergence(m, parameter = "rho"),
    "Items not provided by user. Picking 5 at random."
  )
  expect_error(assess_convergence(m, parameter = "rororo"))

  m <- compute_mallows(potato_visual, nmc = 10, save_aug = FALSE, n_clusters = 2)
  expect_error(assess_convergence(m, parameter = "Rtilde"))

  m <- compute_mallows(preferences = beach_preferences, nmc = 10, error_model = "bernoulli")
  expect_s3_class(assess_convergence(m, parameter = "theta"), "ggplot")

  # This one should not give message
  testdat <- matrix(
    c(
      1, 2, 3,
      2, 1, 3,
      1, 2, 3
    ),
    nrow = 3, ncol = 3, byrow = TRUE
  )

  m <- compute_mallows(testdat, nmc = 5)
  expect_equal(
    assess_convergence(m, parameter = "rho"),
    assess_convergence(m, parameter = "rho", items = 1:3)
  )
  m <- compute_mallows(testdat, nmc = 5, n_clusters = 2, save_aug = TRUE)
  expect_equal(
    assess_convergence(m, parameter = "Rtilde"),
    assess_convergence(m, parameter = "Rtilde", items = 1:3)
  )
})


test_that("assess_convergence works with mixtures", {
  m <- compute_mallows_mixtures(n_clusters = 3, rankings = potato_visual, nmc = 10)
  class(m) <- NULL
  expect_error(assess_convergence(m))
  expect_error(assess_convergence(m, parameter = "aaa"))
  expect_error(assess_convergence(m, parameter = "Rtilde"))

  m <- compute_mallows_mixtures(
    n_clusters = 3, rankings = potato_visual, nmc = 10,
    save_aug = TRUE
  )
  expect_s3_class(assess_convergence(m, parameter = "alpha"), "ggplot")
  expect_error(assess_convergence(m, parameter = "Rtilde"))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = 1:100))
  expect_error(assess_convergence(m, parameter = "Rtilde", assessors = c("a", "b")))

  expect_s3_class(assess_convergence(m), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "rho"), "ggplot")
  expect_s3_class(assess_convergence(m, parameter = "cluster_probs"), "ggplot")

  m <- compute_mallows(potato_visual, nmc = 10, n_clusters = 2, save_clus = TRUE)
  plt <- assess_convergence(m, parameter = "cluster_probs")
  expect_s3_class(plt, "ggplot")
  pdf(NULL)
  print(plt)
  dev.off()
})
