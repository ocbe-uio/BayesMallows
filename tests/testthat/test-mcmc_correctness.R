test_that("compute_mallows is correct for complete data", {
  expectations <- data.frame(
    metric = c("footrule", "spearman", "kendall"),
    mean = c(10.861, 2.37, 16.46),
    sd = c(0.735428, 0.24803, 1.2266)
  )

  for (i in seq_len(nrow(expectations))) {
    set.seed(1)
    mod_bmm <- compute_mallows(
      data = setup_rank_data(potato_visual),
      model_options = set_model_options(metric = expectations$metric[[i]]),
      compute_options = set_compute_options(nmc = 50000, burnin = 1000)
    )

    expect_equal(
      mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$mean[[i]],
      tolerance = .01
    )

    expect_equal(
      sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$sd[[i]],
      tolerance = .1
    )
  }
})

test_that("compute_mallows is correct for pairwise preferences", {
  expectations <- data.frame(
    metric = c("footrule", "kendall"),
    mean = c(4.844578, 6.622525),
    sd = c(0.2916162, 0.4074065)
  )

  for (i in seq_len(nrow(expectations))) {
    set.seed(123)

    mod_bmm <- compute_mallows(
      data = setup_rank_data(preferences = beach_preferences),
      model_options = set_model_options(metric = expectations$metric[[i]]),
      compute_options = set_compute_options(
        nmc = 10000, burnin = 1000, alpha_prop_sd = .1
      )
    )

    expect_equal(
      mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$mean[[i]],
      tolerance = .05
    )

    expect_equal(
      sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$sd[[i]],
      tolerance = .1
    )
  }
})

test_that("augmented rankings obey transitive closure", {
  beach_data <- setup_rank_data(
    preferences = beach_preferences
  )

  model_fit <- compute_mallows(
    data = beach_data,
    compute_options = set_compute_options(save_aug = TRUE)
  )

  expect_true(all(
    subset(model_fit$augmented_data, assessor == 1 & item == "Item 4")$value >
      subset(model_fit$augmented_data, assessor == 1 & item == "Item 1")$value
  ))

  expect_true(all(
    subset(model_fit$augmented_data, assessor == 5 & item == "Item 15")$value >
      subset(model_fit$augmented_data, assessor == 5 & item == "Item 9")$value
  ))

  expect_true(all(
    subset(model_fit$augmented_data, assessor == 15 & item == "Item 7")$value >
      subset(model_fit$augmented_data, assessor == 15 & item == "Item 15")$value
  ))
})

test_that("compute_mallows is correct for top-k ranks", {
  dat <- potato_visual
  dat[dat > 10] <- NA

  expectations <- data.frame(
    metric = c("footrule", "spearman", "kendall"),
    mean = c(9.629511, 1.978859, 14.07414),
    sd = c(0.8694795, 0.3140762, 1.385546)
  )

  for (i in seq_len(nrow(expectations))) {
    set.seed(123)
    mod_bmm <- compute_mallows(
      data = setup_rank_data(dat),
      model_options = set_model_options(metric = expectations$metric[[i]]),
      compute_options = set_compute_options(nmc = 2000, burnin = 500)
    )

    expect_equal(
      mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$mean[[i]],
      tolerance = 1e-4
    )

    expect_equal(
      sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
      expectations$sd[[i]],
      tolerance = 1e-4
    )
  }
})

test_that("compute_mallows is correct with clustering", {
  set.seed(1)
  mixture_model <- compute_mallows(
    data = setup_rank_data(cluster_data),
    model_options = set_model_options(n_clusters = 3),
    compute_options = set_compute_options(nmc = 5000)
  )

  aggdat <- aggregate(value ~ cluster, data = mixture_model$alpha, FUN = mean)
  expect_gt(max(aggdat$value) - min(aggdat$value), 1)

  aggdat <- aggregate(value ~ cluster,
    data = mixture_model$cluster_probs,
    FUN = mean
  )
  expect_true(all(aggdat$value < .5))

  skip_on_cran()
  set.seed(1)
  cl <- parallel::makeCluster(2)
  models <- compute_mallows_mixtures(
    n_clusters = c(1, 5, 10),
    data = setup_rank_data(sushi_rankings),
    compute_options = set_compute_options(nmc = 1000, include_wcd = TRUE),
    cl = cl
  )
  parallel::stopCluster(cl)

  wcd_means <- vapply(models, function(x) {
    mean(x$within_cluster_distance$value[
      x$within_cluster_distance$iteration > 100
    ])
  }, 1)

  expect_equal(
    wcd_means, sort(wcd_means, decreasing = TRUE)
  )

  mixture_model <- compute_mallows(
    data = setup_rank_data(rankings = sushi_rankings),
    compute_options = set_compute_options(nmc = 1000, burnin = 100),
    model_options = set_model_options(n_clusters = 5)
  )

  dat <- mixture_model$alpha[mixture_model$alpha$iteration > 100, ]
  aggdat <- aggregate(value ~ cluster, data = dat, FUN = mean)

  expect_gte(
    max(aggdat$value) - min(aggdat$value),
    1
  )
})

test_that("compute_mallows is correct with Bernoulli error", {
  set.seed(33)
  mod <- compute_mallows(
    setup_rank_data(preferences = bernoulli_data),
    compute_options = set_compute_options(nmc = 2000, swap_leap = 2),
    model_options = set_model_options(error_model = "bernoulli")
  )

  expect_equal(
    mean(mod$alpha$value[mod$alpha$iteration > 1500]),
    9.518804,
    tolerance = 1e-4
  )

  expect_equal(
    sd(mod$alpha$value[mod$alpha$iteration > 1500]),
    0.7038427,
    tolerance = 1e-4
  )

  expect_equal(
    mean(mod$theta$value[mod$theta$iteration > 1500]),
    0.1023917,
    tolerance = 1e-4
  )

  expect_equal(
    sd(mod$theta$value[mod$theta$iteration > 1500]),
    0.004577318,
    tolerance = 1e-4
  )
})

test_that("compute_mallows is correct with Bernoulli error and partial data", {
  set.seed(33)
  mod <- compute_mallows(
    setup_rank_data(
      preferences = bernoulli_data[sample(nrow(bernoulli_data), 500), ]
    ),
    compute_options = set_compute_options(nmc = 10000, swap_leap = 3),
    model_options = set_model_options(error_model = "bernoulli")
  )

  expect_equal(
    mean(mod$alpha$value[mod$alpha$iteration > 8000]),
    9.260634,
    tolerance = 1e-4
  )

  expect_equal(
    sd(mod$alpha$value[mod$alpha$iteration > 8000]),
    0.8246462,
    tolerance = 1e-4
  )

  expect_equal(
    mean(mod$theta$value[mod$theta$iteration > 8000]),
    0.06693134,
    tolerance = 1e-4
  )

  expect_equal(
    sd(mod$theta$value[mod$theta$iteration > 8000]),
    0.009311273,
    tolerance = 1e-4
  )
})

test_that("swap propsal works for modal ranking", {
  # They're not supposed to be identical, but should be fairly close
  # for this easy problem.
  set.seed(3)
  mod1 <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10000, burnin = 2000)
  )

  mod2 <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(
      nmc = 10000, burnin = 2000,
      rho_proposal = "swap"
    )
  )

  expect_equal(
    mean(mod1$alpha$value[mod1$alpha$iteration > 2000]),
    mean(mod2$alpha$value[mod1$alpha$iteration > 2000]),
    tolerance = .01
  )

  expect_equal(
    sd(mod1$alpha$value[mod1$alpha$iteration > 2000]),
    sd(mod2$alpha$value[mod1$alpha$iteration > 2000]),
    tolerance = 1
  )

  expect_equal(
    compute_consensus(mod1)$item,
    compute_consensus(mod2)$item
  )
})
