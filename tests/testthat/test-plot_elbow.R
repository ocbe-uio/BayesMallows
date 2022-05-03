context("Testing plot_elbow")

test_that("plot_elbow fails when it should", {
  test <- compute_mallows_mixtures(
    n_clusters = 1:3, rankings = potato_visual,
    nmc = 20, include_wcd = TRUE
  )
  expect_s3_class(plot_elbow(test, burnin = 16), "ggplot")
  expect_error(plot_elbow(test, burnin = 300))
  expect_error(plot_elbow(test, burnin = 200))
  expect_error(plot_elbow(test))
  test <- compute_mallows_mixtures(
    n_clusters = 1:3, rankings = potato_visual,
    nmc = 20, include_wcd = FALSE
  )
  expect_error(plot_elbow(test, burnin = 100))
})
