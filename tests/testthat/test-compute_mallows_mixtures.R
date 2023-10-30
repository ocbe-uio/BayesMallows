test_that("compute_mallows_mixtures works", {
  set.seed(1234)
  n_clusters <- c(1, 4, 6)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    compute_options = set_compute_options(nmc = 20L),
    rankings = sushi_rankings[1:100, ]
  )

  expect_equal(
    round(models[[1]]$alpha$value, 10),
    c(
      1, 0.9613263178, 0.9613263178, 0.8831344651, 0.9125035114,
      0.9125035114, 0.8865988675, 0.8865988675, 0.8865988675, 0.8865988675,
      0.8205473633, 0.8205473633, 0.6941503586, 0.767450906, 0.7817768491,
      0.882383532, 0.9845249867, 0.9283154749, 1.1219871023, 1.0797258265
    )
  )

  expect_equal(
    models[[2]]$rho$value[12:18],
    c(6, 2, 5, 10, 3, 8, 1)
  )

  expect_equal(
    models[[3]]$rho_acceptance,
    c(0.75, 0.8, 0.95, 0.45, 0.9, 0.45)
  )


  set.seed(123)
  mixture_model <- compute_mallows(
    rankings = sushi_rankings[1:100, ], n_clusters = 5,
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10)
  )

  expect_equal(
    mixture_model$within_cluster_distance$value,
    c(
      684, 740, 716, 770, 506, 726, 574, 434, 420, 856, 678, 502,
      388, 370, 996, 746, 366, 462, 390, 988, 1082, 250, 248, 82, 1170,
      778, 194, 206, 106, 1444, 492, 60, 270, 342, 1580, 618, 70, 106,
      432, 1510, 716, 92, 272, 206, 1402, 728, 156, 196, 170, 1532
    )
  )

  # check that it runs in parallel
  cl <- parallel::makeCluster(1)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    compute_options = set_compute_options(nmc = 20L),
    rankings = sushi_rankings[1:100, ],
    cl = cl
  )
  parallel::stopCluster(cl)
  expect_s3_class(models, "BayesMallowsMixtures")

  # check that psi argument is being used
  set.seed(123)
  mixture_model1 <- compute_mallows(
    rankings = sushi_rankings[1:100, ], n_clusters = 5,
    psi = 100,
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10)
  )

  set.seed(123)
  mixture_model2 <- compute_mallows(
    rankings = sushi_rankings[1:100, ], n_clusters = 5,
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10),
    psi = .1
  )

  expect_lt(
    max(mixture_model1$cluster_probs$value),
    max(mixture_model$cluster_probs$value)
  )

  expect_gt(
    max(mixture_model2$cluster_probs$value),
    max(mixture_model$cluster_probs$value)
  )
})
