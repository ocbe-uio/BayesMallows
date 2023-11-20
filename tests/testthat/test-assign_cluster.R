test_that("assign_cluster works", {
  m <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10)
  )

  expect_error(assign_cluster(m))

  set.seed(123)
  m <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10)
  )

  asc <- assign_cluster(m, burnin = 5)
  expect_equal(unique(asc$cluster), "Cluster 1")
  expect_equal(unique(asc$probability), 1)
  expect_equal(unique(asc$map_cluster), "Cluster 1")

  set.seed(123)
  m <- compute_mallows(
    data = setup_rank_data(potato_visual),
    model = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 10)
  )

  asc <- assign_cluster(m, burnin = 7)
  asc <- asc[order(as.integer(asc$assessor)), ]

})
