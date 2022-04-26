test_that("assign_cluster works", {

  m <- compute_mallows(potato_visual, nmc = 10)

  expect_error(assign_cluster(m))
  expect_error(assign_cluster(m, burnin = 1))

  set.seed(123)
  m <- compute_mallows(potato_visual, nmc = 10, save_clus = TRUE)

  asc <- assign_cluster(m, burnin = 5)
  expect_equal(unique(asc$cluster), "Cluster 1")
  expect_equal(unique(asc$probability), 1)
  expect_equal(unique(asc$map_cluster), "Cluster 1")

  set.seed(123)
  m <- compute_mallows(potato_visual, nmc = 10, save_clus = TRUE, n_clusters = 3)

  asc <- assign_cluster(m, burnin = 7)
  expect_equal(
    round(asc$probability, 4),
    c(0.6667, 0.3333, 1, 0.6667, 0.3333, 0.3333, 0.3333, 0.3333,
      0.6667, 0.3333, 1, 0.6667, 0.3333, 1, 0.3333, 0.3333, 0.3333,
      0.3333, 0.3333, 0.3333, 1, 1)
  )

  expect_equal(
    asc$map_cluster,
    c("Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 2")
  )
  expect_equal(
    asc$cluster,
    c("Cluster 1", "Cluster 2", "Cluster 1", "Cluster 1", "Cluster 2",
      "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 3",
      "Cluster 1", "Cluster 1", "Cluster 2", "Cluster 1", "Cluster 1",
      "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2", "Cluster 3",
      "Cluster 1", "Cluster 2")
  )

  asc <- assign_cluster(m, burnin = 9, expand = TRUE)
  expect_equal(
    round(asc$probability[order(as.integer(asc$assessor))], 4),
    c(0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0,
      0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0)
    )

  expect_equal(
    asc$cluster[order(as.integer(asc$assessor))],
    c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2",
      "Cluster 3", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1",
      "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2", "Cluster 3",
      "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2",
      "Cluster 3", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1",
      "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2", "Cluster 3",
      "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 1", "Cluster 2",
      "Cluster 3")
  )

  expect_equal(
    asc$map_cluster[order(as.integer(asc$assessor))],
    c("Cluster 2", "Cluster 2", "Cluster 2", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 2", "Cluster 2", "Cluster 2", "Cluster 2",
      "Cluster 2", "Cluster 2", "Cluster 3", "Cluster 3", "Cluster 3",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1",
      "Cluster 1", "Cluster 1", "Cluster 3", "Cluster 3", "Cluster 3",
      "Cluster 1", "Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2",
      "Cluster 2")
  )

  asc <- assign_cluster(m, burnin = 8, soft = FALSE)
  expect_equal(dim(asc), c(12L, 3L))
  expect_equal(asc$probability,
               c(0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 1, 0.5, 0.5, 1, 1))

})
