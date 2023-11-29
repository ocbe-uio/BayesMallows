test_that("compute_mallows is correct for complete data", {
  expectations <- data.frame(
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
    mean = c(10.85, 1.988, 20.11, 29.53, 16.43, 35.03),
    sd = c(0.77, 0.197, 4.31, 3.10, 1.287, 4.149)
  )

  for (i in seq_len(nrow(expectations))) {
    set.seed(123)
    mod_bmm <- compute_mallows(
      data = setup_rank_data(potato_visual),
      model_options = set_model_options(metric = expectations$metric[[i]]),
      compute_options = set_compute_options(nmc = 10000, burnin = 1000)
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


test_that("compute_mallows is correct for pairwise preferences", {
  expectations <- data.frame(
    metric = c("footrule", "kendall"),
    mean = c(4.825, 6.666),
    sd = c(0.286, 0.418)
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
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
    mean = c(10.033, 1.686, 35.28, 44.637, 14.70, 52.09),
    sd = c(0.781, 0.6463, 5.91, 3.44, 1.38, 5.19)
  )

  for (i in seq_len(nrow(expectations))) {
    set.seed(123)
    mod_bmm <- compute_mallows(
      data = setup_rank_data(dat),
      model_options = set_model_options(metric = expectations$metric[[i]]),
      compute_options = set_compute_options(nmc = 200000, burnin = 100000)
    )

    expect_equal(
      mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 100000]),
      expectations$mean[[i]],
      tolerance = .05
    )

    expect_equal(
      sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 10000]),
      expectations$sd[[i]],
      tolerance = .2
    )
  }
})
