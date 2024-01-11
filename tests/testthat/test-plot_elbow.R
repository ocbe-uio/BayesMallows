test_that("plot_elbow works", {
  set.seed(1)
  n_clusters <- seq(from = 1, to = 5)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(cluster_data),
    compute_options = set_compute_options(nmc = 2, include_wcd = TRUE)
  )

  expect_error(plot_elbow(models), "burnin not provided")

  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(cluster_data),
    compute_options =
      set_compute_options(nmc = 3, burnin = 0, include_wcd = FALSE)
  )

  expect_error(
    plot_elbow(models),
    "To get an elbow plot, set include_wcd=TRUE"
  )

  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(cluster_data),
    compute_options =
      set_compute_options(nmc = 10, burnin = 1, include_wcd = TRUE)
  )

  p <- plot_elbow(models)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Within-cluster sum of distances")
  expect_equal(p$labels$x, "Number of clusters")
  expect_equal(dim(p$data), c(45, 3))

  mod <- compute_mallows(
    data = setup_rank_data(cluster_data),
    model_options = set_model_options(metric = "cayley", n_clusters = 3),
    compute_options = set_compute_options(nmc = 100, include_wcd = TRUE)
  )

  expect_error(plot_elbow(mod), "burnin not provided")
  mod$burnin <- 200
  expect_error(
    plot_elbow(mod),
    "burnin must be strictly smaller than the number of MCMC samples"
  )
  mod$burnin <- 0
  expect_s3_class(plot_elbow(mod), "ggplot")
})
