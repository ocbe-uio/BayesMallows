cl <- parallel::makeCluster(2)
set.seed(1234)

test_that("miscellaneous input validation", {
  namat <- potato_visual
  namat[c(1, 2, 3), c(7, 9)] <- NA_real_
  expect_error(compute_mallows(rankings = namat, na_action = "fail", cl = cl))
  expect_output(
    compute_mallows(rankings = namat,
                    model = set_model_options(na_action = "omit"),
                    compute_options = set_compute_options(nmc = 2)),
    "Omitting 9 rows from rankings due to NA values"
  )
  expect_s3_class(
    compute_mallows(
      rankings = namat,
      model = set_model_options(na_action = "augment"),
      compute_options = set_compute_options(nmc = 3),
      cl = cl),
    "BayesMallows")
  expect_s3_class(
    compute_mallows(rankings = namat,
                    compute_options = set_compute_options(nmc = 3),
                    cl = cl),
    "BayesMallows")
  expect_error(
    compute_mallows(
      compute_options = set_compute_options(nmc = 1000, alpha_prop_sd = 1),
      cl = cl))
  expect_error(
    compute_mallows(rankings = potato_visual,
                    compute_options = set_compute_options(nmc = 100, alpha_jump = 102), cl = cl))
  expect_error(
    compute_mallows(rankings = potato_visual, lambda = 0, cl = cl))
  expect_error(
    compute_mallows(rankings = potato_visual, lambda = -10, cl = cl))
  expect_error(compute_mallows(rankings = potato_visual,
                               compute_options = set_compute_options(nmc = 100, rho_thinning = 200), cl = cl))
  expect_error(compute_mallows(
    rankings = potato_visual,
    compute_options = set_compute_options(nmc = 100, aug_thinning = 200), cl = cl))
  expect_error(
    compute_mallows(rankings = potato_visual,
                    compute_options = set_compute_options(nmc = -100), cl = cl))
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
  expect_error(invisible(capture.output(
    compute_mallows(preferences = preferences,
                    compute_options = set_compute_options(nmc = 10), cl = cl))))
  expect_s3_class(
    compute_mallows(
      preferences = preferences,
      model = set_model_options(error_model = "bernoulli"),
      compute_options = set_compute_options(nmc = 10)),
    "BayesMallows"
  )
})

test_that("compute_mallows with single missing value works", {
  dd <- potato_visual
  dd[1, 1] <- NA
  dd[2, 3] <- NA
  m <- compute_mallows(dd, compute_options = set_compute_options(nmc = 4),
                       cl = cl)
  expect_gt(mean(m$alpha$value), 0)
})

test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(rankings = mat,
                       compute_options = set_compute_options(nmc = 30), cl = cl)
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")
})


test_that("compute_mallows runs with the right distances", {
  for (metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")) {
    expect_s3_class(
      compute_mallows(potato_visual,
                      model = set_model_options(metric = metric),
                      compute_options = set_compute_options(nmc = 3), cl = cl), "BayesMallows")
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

  expect_s3_class(
    compute_mallows(preferences = m,
                    compute_options = set_compute_options(nmc = 20), cl = cl),
    "BayesMallows")
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
    rankings = potato_visual, compute_options = set_compute_options(nmc = 10),
    obs_freq = rep(1, nrow(potato_visual)), seed = 1L, cl = cl
  )
  m2 <- compute_mallows(
    rankings = potato_visual, compute_options = set_compute_options(nmc = 10),
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
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    seed = 3344L, cl = cl
  )

  expect_equal(
    model_fit_obs_freq$rho$value,
    c(14, 5, 4, 6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 5, 3,
      6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 4, 14, 5, 4, 6, 11, 8,
      12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 5, 3, 6, 11, 8, 12, 10, 1,
      7, 15, 13, 9, 2, 4, 14, 5, 3, 6, 12, 8, 11, 10, 1, 7, 15, 13,
      9, 2, 4, 14, 5, 4, 6, 12, 8, 11, 10, 1, 7, 15, 13, 9, 2, 3, 14,
      5, 4, 6, 12, 8, 11, 10, 1, 7, 15, 13, 9, 2, 3, 14, 5, 4, 6, 12,
      9, 11, 10, 1, 7, 15, 13, 8, 2, 3, 14, 5, 4, 6, 12, 9, 10, 11,
      1, 7, 15, 13, 8, 2, 3, 14, 5, 4, 6, 12, 9, 10, 11, 1, 7, 15,
      13, 8, 3, 2, 12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14,
      12, 13, 9, 3, 5, 15, 10, 1, 4, 8, 7, 11, 2, 6, 14, 12, 13, 10,
      3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14, 12, 13, 10, 3, 5, 15,
      9, 1, 4, 8, 6, 11, 2, 7, 14, 12, 13, 10, 3, 6, 15, 9, 1, 4, 8,
      5, 11, 2, 7, 14, 12, 13, 10, 3, 6, 15, 9, 1, 4, 8, 5, 11, 2,
      7, 14, 12, 14, 10, 3, 6, 15, 9, 1, 4, 8, 5, 11, 2, 7, 13, 12,
      14, 10, 3, 6, 15, 9, 1, 4, 8, 5, 11, 2, 7, 13, 12, 14, 9, 3,
      6, 15, 10, 1, 4, 8, 5, 11, 2, 7, 13, 12, 14, 9, 2, 6, 15, 10,
      1, 4, 8, 5, 11, 3, 7, 13)
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
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    seed = 3344L, cl = cl
  )

  expect_equal(
    model_fit_rep$rho$value,
    c(14, 5, 4, 6, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 6, 4,
      5, 11, 8, 12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 6, 5, 4, 11, 8,
      12, 10, 1, 7, 15, 13, 9, 2, 3, 14, 6, 5, 4, 11, 7, 12, 10, 1,
      8, 15, 13, 9, 2, 3, 14, 6, 5, 4, 11, 7, 12, 10, 1, 8, 15, 13,
      9, 3, 2, 14, 6, 5, 4, 11, 7, 12, 10, 1, 8, 15, 13, 9, 3, 2, 14,
      6, 5, 4, 11, 7, 13, 10, 1, 8, 15, 12, 9, 3, 2, 15, 6, 5, 4, 11,
      7, 13, 10, 1, 8, 14, 12, 9, 3, 2, 15, 6, 5, 4, 11, 7, 13, 10,
      1, 9, 14, 12, 8, 3, 2, 15, 6, 5, 4, 11, 7, 13, 10, 1, 9, 14,
      12, 8, 3, 2, 12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14,
      12, 13, 10, 3, 5, 15, 9, 1, 4, 8, 7, 11, 2, 6, 14, 12, 13, 11,
      3, 5, 15, 9, 1, 4, 8, 7, 10, 2, 6, 14, 12, 13, 11, 3, 5, 15,
      9, 1, 4, 8, 7, 10, 2, 6, 14, 12, 13, 11, 3, 4, 15, 9, 1, 5, 8,
      7, 10, 2, 6, 14, 12, 13, 11, 3, 4, 15, 10, 1, 5, 8, 7, 9, 2,
      6, 14, 12, 14, 11, 3, 4, 15, 10, 1, 5, 8, 7, 9, 2, 6, 13, 12,
      14, 10, 3, 4, 15, 11, 1, 5, 8, 7, 9, 2, 6, 13, 13, 14, 10, 3,
      4, 15, 11, 1, 5, 8, 7, 9, 2, 6, 12, 13, 14, 10, 3, 5, 15, 11,
      1, 4, 8, 7, 9, 2, 6, 12)
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
  mm <- compute_mallows(potato_visual,
                        compute_options = set_compute_options(nmc = 5),
                        init = set_initial_values(alpha_init = c(3, 4)),
                        cl = cl)
  expect_s3_class(mm, "BayesMallows")
})

parallel::stopCluster(cl)
