test_that("compute_mallows_mixtures works", {
  n_clusters <- c(1, 4, 6)
  dat <- setup_rank_data(rankings = sushi_rankings[1:100, ])

  set.seed(1234)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    data = dat,
    compute_options = set_compute_options(nmc = 20L)
  )

  expect_equal(
    round(models[[1]]$alpha$value, 10),
    c(
      1, 0.9613263178, 0.9613263178, 0.8831344651, 0.8831344651,
      0.8831344651, 0.8580635656, 0.8580635656, 0.8580635656, 0.8580635656,
      0.7941379378, 0.7941379378, 0.6718090374, 0.6718090374, 0.6843496416,
      0.6843496416, 0.6843496416, 0.645278049, 0.645278049, 0.6209727129
    )
  )

  expect_equal(
    models[[2]]$rho$value[12:18],
    c(6, 2, 5, 10, 3, 8, 1)
  )

  set.seed(123)
  mixture_model <- compute_mallows(
    data = dat,
    model_options = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10)
  )


  # check that it runs in parallel
  cl <- parallel::makeCluster(1)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    data = dat,
    compute_options = set_compute_options(nmc = 20L),
    cl = cl
  )
  parallel::stopCluster(cl)
  expect_s3_class(models, "BayesMallowsMixtures")

  # check that psi argument is being used
  set.seed(123)
  mixture_model1 <- compute_mallows(
    data = dat,
    model_options = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 1000),
    priors = set_priors(psi = 100)
  )

  set.seed(123)
  mixture_model2 <- compute_mallows(
    data = dat,
    model_options = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 1000),
    priors = set_priors(psi = 1)
  )

  expect_lt(
    max(mixture_model1$cluster_probs$value),
    max(mixture_model2$cluster_probs$value)
  )
})
