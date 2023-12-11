test_that("heat_plot works", {
  set.seed(1)
  model_fit <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 2)
  )
  expect_s3_class(heat_plot(model_fit, burnin = 0, type = "MAP"), "ggplot")

  set.seed(1)
  model_fit <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 20)
  )

  expect_error(heat_plot(model_fit), "Please specify the burnin.")
  expect_error(heat_plot(model_fit, burnin = 300), "burnin must be <= nmc")

  expect_s3_class(heat_plot(model_fit, burnin = 1), "ggplot")
  model_fit$burnin <- 1
  expect_s3_class(heat_plot(model_fit, type = "MAP"), "ggplot")

  model_fit <- compute_mallows(
    setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 20, burnin = 10)
  )
  expect_error(
    heat_plot(model_fit),
    "heat_plot only works for a single cluster"
  )
})
