test_that("assign_cluster fails properly", {
  mod <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10)
    )

  expect_error(assign_cluster(mod), "Please specify the burnin.")
  mod$burnin <- 11
  expect_error(assign_cluster(mod), "burnin < model_fit")

})

test_that("assign_cluster works", {
  set.seed(123)
  mod <- compute_mallows(
    setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 300, burnin = 50)
    )

  a1 <- assign_cluster(mod, soft = FALSE, expand = FALSE)
  expect_equal(dim(a1), c(60, 3))
  agg1 <- aggregate(assessor ~ map_cluster, a1, length)
  expect_equal(agg1$assessor, c(17, 24, 19))

  a2 <- assign_cluster(mod, soft = TRUE, expand = FALSE)
  expect_equal(ncol(a2), 4)
  agg2 <- aggregate(probability ~ assessor, a2, sum)
  expect_equal(mean(agg2$probability), 1)

  expect_equal(
    dim(assign_cluster(mod, soft = FALSE, expand = TRUE)),
    c(60, 3))

  a3 <- assign_cluster(mod, soft = TRUE, expand = TRUE)
  agg3 <- aggregate(probability ~ assessor, a3, sum)
  expect_equal(mean(agg2$probability), 1)

  mod <- compute_mallows(
    setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 2, burnin = 1)
  )

  expect_equal(dim(assign_cluster(mod)), c(60, 4))
  expect_equal(dim(assign_cluster(mod, expand = TRUE)), c(180, 4))

})
