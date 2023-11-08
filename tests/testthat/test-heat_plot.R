model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 200), seed = 1
)


test_that("heat_plot fails when it should", {
  expect_error(heat_plot(model_fit))
  expect_error(heat_plot(model_fit, burnin = 1e8))
  expect_error(heat_plot(compute_mallows(potato_visual, n_clusters = 2), burnin = 10))
})

model_fit$burnin <- 10

expect_s3_class(heat_plot(model_fit), "ggplot")
expect_s3_class(heat_plot(model_fit, type = "MAP"), "ggplot")
