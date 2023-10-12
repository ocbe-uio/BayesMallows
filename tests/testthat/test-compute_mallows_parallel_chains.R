context("Testing compute_mallows with parallel chains")
cl <- parallel::makeCluster(2)
set.seed(1234)

test_that("miscellaneous input validation", {
  namat <- potato_visual
  namat[c(1, 2, 3), c(7, 9)] <- NA_real_
  expect_error(compute_mallows(rankings = namat, na_action = "fail", cl = cl))
  expect_output(
    compute_mallows(rankings = namat, nmc = 2, na_action = "omit"),
    "Omitting 9 rows from rankings due to NA values"
  )
  expect_s3_class(compute_mallows(rankings = namat, na_action = "augment", nmc = 3, cl = cl), "BayesMallows")
  expect_s3_class(compute_mallows(rankings = namat, nmc = 3, cl = cl), "BayesMallows")
  expect_error(compute_mallows(nmc = 1000, alpha_prop_sd = 1, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, alpha_jump = 102, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, lambda = 0, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, lambda = -10, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, rho_thinning = 200, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, aug_thinning = 200, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual, nmc = -100, cl = cl))
})

test_that("rho_init is properly validated", {
  m <- potato_visual
  expect_error(compute_mallows(rankings = m, rho_init = 1:(ncol(m) - 1), cl = cl))
  expect_error(compute_mallows(rankings = m, rho_init = c(potato_true_ranking[-1], 22), cl = cl))
  expect_error(compute_mallows(rankings = m, rho_init = c(NA_real_, 2:ncol(m)), cl = cl))
  expect_error(compute_mallows(rankings = m, obs_freq = -1, cl = cl))
  expect_error(compute_mallows(rankings = m, obs_freq = 1, cl = cl))
  expect_error(compute_mallows(rankings = m, obs_freq = 1:11, cl = cl))
})

test_that("compute_mallows discovers inconsistent rankings", {
  expect_error(compute_mallows(
    rankings = matrix(c(
      1, 2, -3,
      1, 2, 3
    ), nrow = 2, byrow = TRUE), cl = cl
  ))
  expect_error(compute_mallows(
    rankings = matrix(c(
      1, 2, 3,
      1, 2, 2
    ), nrow = 2, byrow = TRUE), cl = cl
  ))
})


test_that("compute_mallows error model works", {
  preferences <- data.frame(
    assessor = c(1, 1, 2, 2),
    bottom_item = c(1, 2, 1, 2),
    top_item = c(2, 1, 2, 3)
  )
  expect_error(invisible(capture.output(compute_mallows(preferences = preferences, nmc = 10, cl = cl))))
  expect_s3_class(
    compute_mallows(preferences = preferences, error_model = "bernoulli", nmc = 10),
    "BayesMallows"
  )
})

test_that("compute_mallows with single missing value works", {
  dd <- potato_visual
  dd[1, 1] <- NA
  dd[2, 3] <- NA
  m <- compute_mallows(dd, nmc = 4, cl = cl)
  expect_gt(mean(m$alpha$value), 0)
})

test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(rankings = mat, nmc = 30, cl = cl)
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")
})


test_that("compute_mallows runs with the right distances", {
  for (metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")) {
    expect_s3_class(compute_mallows(potato_visual, metric = metric, nmc = 3, cl = cl), "BayesMallows")
  }
})

test_that("compute_mallows handles integer preferences", {
  set.seed(123)
  m <- subset(
    beach_preferences,
    top_item %in% c(1, 2, 3) | bottom_item %in% c(1, 2, 3)
  )
  m[sample(nrow(m), 20), , ]
  for (col in names(m)) {
    eval(parse(text = paste("m$", col, "<- as.integer(m$", col, ")")))
  }

  expect_s3_class(compute_mallows(preferences = m, nmc = 20, cl = cl), "BayesMallows")
})

test_that("compute_mallows handles data with lots of missings", {
  R_partial2 <- structure(c(
    NA, NA, NA, NA, NA, NA, 9, NA, NA, 7, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, 7, 8, 10, NA, NA, 9, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6, NA, 5, 2,
    6, 5, 6, 6, 5, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, 3, 4, NA, 3, 3, 3, 3, 4, 5, 3, 3, 3, NA, 3, 3, 4, NA,
    7, 8, 3, 3, 10, 5, 4, NA, NA, NA, 8, NA, NA, NA, NA, NA, 11,
    NA, NA, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 4, 2,
    2, 2, 4, 2, 2, 2, NA, NA, 4, 7, 5, 4, 6, 7, 2, 6, 6, 7, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9,
    7, 8, NA, 8, 10, 6, NA, 5, NA, 6, 6, 5, 4, 5, NA, 4, 4, 5, NA,
    NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, 9, NA
  ), .Dim = c(12L, 20L))

  m <- compute_mallows(R_partial2, cl = cl)
  expect_s3_class(assess_convergence(m), "gg")
})

test_that("compute_mallows treats obs_freq properly", {
  m1 <- compute_mallows(
    rankings = potato_visual, nmc = 10,
    obs_freq = rep(1, nrow(potato_visual)), seed = 1L, cl = cl
  )
  m2 <- compute_mallows(
    rankings = potato_visual, nmc = 10,
    seed = 1L, cl = cl
  )
  expect_equal(m1, m2)

  # Test with repeated beach preferences
  obs_freq <- c(2, 1, 4)

  beach_small <- subset(beach_preferences, assessor %in% c(1, 2, 3))

  # Next, we create a new hypothetical beach_preferences dataframe where each
  # assessor is replicated 1-4 times

  beach_pref_rep <- do.call(rbind, lapply(split(beach_small, f = seq_len(nrow(beach_small))), function(dd) {
    ret <- merge(
      dd,
      data.frame(new_assessor = seq_len(obs_freq[dd$assessor])),
      all = TRUE
    )
    ret$assessor <- paste(ret$assessor, ret$new_assessor, sep = ",")
    ret$new_assessor <- NULL
    ret
  }))


  # We generate transitive closure for these preferences
  beach_tc_rep <- generate_transitive_closure(beach_pref_rep)

  # We generate the initial rankings for the repeated and the "unrepeated"
  # data
  set.seed(1223)
  beach_tc <- generate_transitive_closure(beach_small)
  beach_rankings <- generate_initial_ranking(beach_tc, n_items = 15)
  beach_rankings_rep <- generate_initial_ranking(beach_tc_rep, n_items = 15)

  model_fit_obs_freq <- compute_mallows(
    rankings = beach_rankings,
    preferences = beach_tc,
    obs_freq = obs_freq,
    save_aug = TRUE,
    nmc = 10, seed = 3344L, cl = cl
  )

  expect_equal(
    model_fit_obs_freq$rho$value,
    c(
      14, 5, 4, 6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 5, 2,
      6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 3, 4, 14, 5, 2, 6, 11, 8,
      12, 10, 1, 7, 15, 13, 9, 4, 3, 14, 5, 2, 6, 11, 8, 12, 10, 1,
      7, 15, 13, 9, 4, 3, 14, 5, 2, 6, 12, 8, 11, 10, 1, 7, 15, 13,
      9, 4, 3, 14, 5, 3, 6, 12, 8, 11, 10, 1, 7, 15, 13, 9, 4, 2, 14,
      5, 4, 6, 12, 8, 11, 10, 2, 7, 15, 13, 9, 1, 3, 14, 5, 4, 6, 12,
      9, 11, 10, 2, 8, 15, 13, 7, 1, 3, 14, 5, 4, 6, 12, 10, 8, 11,
      2, 9, 15, 13, 7, 1, 3, 14, 5, 4, 6, 12, 10, 8, 11, 1, 9, 15,
      13, 7, 3, 2, 12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14,
      12, 13, 7, 3, 5, 15, 10, 1, 4, 9, 8, 11, 2, 6, 14, 12, 13, 8,
      3, 5, 15, 7, 1, 4, 10, 9, 11, 2, 6, 14, 12, 13, 8, 3, 5, 15,
      7, 1, 4, 10, 9, 11, 2, 6, 14, 12, 13, 9, 3, 5, 15, 7, 1, 4, 10,
      8, 11, 2, 6, 14, 12, 13, 9, 3, 5, 15, 7, 1, 4, 10, 8, 11, 2,
      6, 14, 13, 14, 9, 3, 5, 15, 7, 1, 4, 10, 8, 12, 2, 6, 11, 13,
      14, 9, 3, 5, 15, 7, 1, 4, 10, 8, 12, 2, 6, 11, 13, 14, 8, 3,
      5, 15, 9, 1, 4, 10, 7, 12, 2, 6, 11, 13, 14, 8, 2, 5, 15, 9,
      1, 3, 10, 7, 12, 4, 6, 11
    )
  )

  expect_equal(
    model_fit_obs_freq$alpha$value,
    c(
      1, 1.06375201942709, 1.23638446536771, 1.18903533953607, 1.12499350243193,
      0.983202718639466, 0.970847993206225, 1.02691791700127, 1.05732241361015,
      0.848313407171291, 1, 1, 1, 1, 1, 1, 0.910488058825459, 0.818043766566705,
      0.818043766566705, 0.743101499728633
    )
  )

  # Next for the repeated data.
  model_fit_rep <- compute_mallows(
    rankings = beach_rankings_rep,
    preferences = beach_tc_rep,
    save_aug = TRUE,
    nmc = 10, seed = 3344L, cl = cl
  )


  expect_equal(
    model_fit_rep$rho$value,
    c(
      14, 5, 4, 6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 8, 4,
      5, 11, 7, 12, 10, 1, 6, 15, 13, 9, 2, 3, 14, 8, 5, 2, 11, 7,
      12, 10, 1, 6, 15, 13, 9, 3, 4, 14, 8, 5, 2, 11, 6, 12, 10, 1,
      7, 15, 13, 9, 3, 4, 14, 8, 5, 3, 11, 6, 12, 10, 1, 7, 15, 13,
      9, 4, 2, 14, 8, 5, 3, 11, 6, 12, 10, 1, 7, 15, 13, 9, 4, 2, 14,
      8, 5, 3, 12, 6, 13, 10, 1, 7, 15, 11, 9, 4, 2, 15, 8, 5, 3, 12,
      6, 13, 10, 1, 7, 14, 11, 9, 4, 2, 15, 7, 5, 3, 12, 6, 13, 9,
      1, 10, 14, 11, 8, 4, 2, 15, 7, 5, 3, 12, 6, 13, 9, 1, 10, 14,
      11, 8, 4, 2, 12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14,
      12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14, 11, 13, 12,
      3, 5, 15, 9, 1, 4, 8, 7, 10, 2, 6, 14, 11, 13, 12, 3, 5, 15,
      9, 1, 4, 8, 7, 10, 2, 6, 14, 11, 13, 12, 4, 2, 15, 9, 1, 5, 8,
      7, 10, 3, 6, 14, 10, 13, 11, 4, 2, 15, 12, 1, 5, 8, 7, 9, 3,
      6, 14, 10, 14, 11, 4, 2, 15, 12, 1, 5, 8, 7, 9, 3, 6, 13, 11,
      14, 9, 4, 2, 15, 12, 1, 5, 8, 7, 10, 3, 6, 13, 12, 14, 9, 4,
      2, 15, 13, 1, 5, 8, 7, 11, 3, 6, 10, 12, 14, 9, 3, 5, 15, 13,
      1, 4, 8, 7, 11, 2, 6, 10
    )
  )

  expect_equal(
    model_fit_rep$alpha$value,
    c(
      1, 1, 1.19997449202174, 1.04873333072888, 1.15559892469591,
      1.08094257108279, 1.07869198268456, 1.05402584033873, 1.20165043028671,
      1.30638245220182, 1, 0.970031407694195, 0.908995832092313, 0.899961691712196,
      0.922780739528762, 0.922780739528762, 0.878378621847317, 0.878286218390065,
      0.755742998802395, 0.755742998802395
    )
  )
})

test_that("compute_mallows() takes initial alpha", {
  mm <- compute_mallows(potato_visual, nmc = 5, cl = cl, alpha_init = c(3, 4))
  expect_s3_class(mm, "BayesMallows")
})

parallel::stopCluster(cl)
