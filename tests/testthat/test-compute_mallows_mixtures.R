test_that("compute_mallows_mixtures works", {
  n_clusters <- c(1, 4, 6)
  dat <- setup_rank_data(rankings = sushi_rankings[1:100, ])

  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    data = dat,
    compute_options = set_compute_options(nmc = 20L),
    seed = 1234
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
    c(10, 8, 3, 7, 5, 2, 1)
  )

  mixture_model <- compute_mallows(
    data = dat,
    model = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10),
    seed = 123
  )

  expect_equal(
    mixture_model$within_cluster_distance$value,
    c(
      684, 740, 716, 770, 506, 726, 574, 434, 420, 856, 678, 672,
      370, 244, 996, 662, 462, 450, 294, 1070, 916, 218, 288, 84, 1334,
      572, 128, 172, 192, 1682, 548, 58, 258, 270, 1564, 632, 184,
      70, 354, 1506, 652, 152, 174, 266, 1526, 772, 144, 166, 318,
      1402
    )
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
  mixture_model1 <- compute_mallows(
    data = dat,
    model = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10),
    priors = set_priors(psi = 100),
    seed = 123
  )

  mixture_model2 <- compute_mallows(
    data = dat,
    model = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE, nmc = 10),
    priors = set_priors(psi = 1),
    seed = 123
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
