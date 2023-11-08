cl <- parallel::makeCluster(2)
set.seed(1234)

test_that("miscellaneous input validation", {
  namat <- potato_visual
  namat[c(1, 2, 3), c(7, 9)] <- NA_real_
  expect_error(
    setup_rank_data(namat, na_action = "fail"),
    "rankings matrix contains NA values"
  )
  expect_output(
    setup_rank_data(namat, na_action = "omit"),
    "Omitting 9 rows from rankings due to NA values"
  )

  expect_s3_class(
    compute_mallows(
      data = setup_rank_data(rankings = namat, na_action = "augment"),
      compute_options = set_compute_options(nmc = 3),
      cl = cl
    ),
    "BayesMallows"
  )
  expect_s3_class(
    compute_mallows(
      data = setup_rank_data(rankings = namat),
      compute_options = set_compute_options(nmc = 3),
      cl = cl
    ),
    "BayesMallows"
  )
  expect_error(
    compute_mallows(
      compute_options = set_compute_options(nmc = 1000, alpha_prop_sd = 1),
      cl = cl
    ),
    "argument \"data\" is missing, with no default"
  )
  expect_error(
    compute_mallows(setup_rank_data(rankings = potato_visual),
      compute_options = set_compute_options(nmc = 100, alpha_jump = 102), cl = cl
    ),
    "nmc must be strictly larger than alpha_jump"
  )
  expect_error(
    compute_mallows(
      setup_rank_data(rankings = potato_visual),
      priors = set_priors(lambda = 0),
      cl = cl
    ),
    "lambda must be a strictly positive number of length one"
  )
  expect_error(
    compute_mallows(
      setup_rank_data(rankings = potato_visual),
      priors = set_priors(lambda = -10),
      cl = cl
    ),
    "lambda must be a strictly positive number of length one"
  )


  expect_error(
    compute_mallows(
      setup_rank_data(rankings = potato_visual),
      compute_options = set_compute_options(nmc = 100, rho_thinning = 200),
      cl = cl
    ),
    "nmc must be strictly larger than rho_thinning"
  )

  expect_error(
    compute_mallows(
      setup_rank_data(rankings = potato_visual),
      compute_options = set_compute_options(nmc = 100, aug_thinning = 200),
      cl = cl
    ),
    "nmc must be strictly larger than aug_thinning"
  )

  expect_error(
    compute_mallows(
      setup_rank_data(rankings = potato_visual),
      compute_options = set_compute_options(nmc = -100),
      cl = cl
    ),
    "nmc must be a positive integer"
  )
})

test_that("rho_init is properly validated", {
  dat <- setup_rank_data(potato_visual)
  expect_error(compute_mallows(data = dat, rho_init = 1:(ncol(m) - 1), cl = cl))
  expect_error(compute_mallows(data = dat, rho_init = c(potato_true_ranking[-1], 22), cl = cl))
  expect_error(compute_mallows(data = dat, rho_init = c(NA_real_, 2:ncol(m)), cl = cl))
  expect_error(compute_mallows(data = dat, observation_frequency = -1, cl = cl))
  expect_error(compute_mallows(data = dat, observation_frequency = 1, cl = cl))
  expect_error(compute_mallows(data = dat, observation_frequency = 1:11, cl = cl))
})

test_that("compute_mallows discovers inconsistent rankings", {
  expect_error(
    setup_rank_data(
      rankings = matrix(c(
        1, 2, -3,
        1, 2, 3
      ), nrow = 2, byrow = TRUE), cl = cl
    ),
    "invalid permutations provided in rankings matrix"
  )
  expect_error(
    setup_rank_data(
      rankings = matrix(c(
        1, 2, 3,
        1, 2, 2
      ), nrow = 2, byrow = TRUE), cl = cl
    ),
    "invalid permutations provided in rankings matrix"
  )
})

test_that("compute_mallows error model works", {
  preferences <- data.frame(
    assessor = c(1, 1, 2, 2),
    bottom_item = c(1, 2, 1, 2),
    top_item = c(2, 1, 2, 3)
  )
  dat <- setup_rank_data(preferences = preferences)
  expect_error(
    compute_mallows(data = dat, cl = cl),
    "Intransitive pairwise comparisons. Please specify an error model."
  )
  expect_s3_class(
    compute_mallows(
      data = dat,
      model = set_model_options(error_model = "bernoulli"),
      compute_options = set_compute_options(nmc = 10),
      cl = cl
    ),
    "BayesMallows"
  )
})

test_that("compute_mallows with single missing value works", {
  dd <- potato_visual
  dd[1, 1] <- NA
  dd[2, 3] <- NA

  m <- compute_mallows(
    setup_rank_data(dd),
    compute_options = set_compute_options(nmc = 4),
    cl = cl
  )
  expect_gt(mean(m$alpha$value), 0)
})

test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(
    setup_rank_data(mat),
    compute_options = set_compute_options(nmc = 30), cl = cl
  )
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")
})


test_that("compute_mallows runs with the right distances", {
  dat <- setup_rank_data(potato_visual)
  for (metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")) {
    expect_s3_class(
      compute_mallows(
        data = dat,
        model = set_model_options(metric = metric),
        compute_options = set_compute_options(nmc = 3), cl = cl
      ), "BayesMallows"
    )
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

  dat <- setup_rank_data(preferences = m)
  expect_s3_class(
    compute_mallows(setup_rank_data(preferences = m),
      compute_options = set_compute_options(nmc = 20), cl = cl
    ),
    "BayesMallows"
  )
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

  m <- compute_mallows(setup_rank_data(R_partial2), cl = cl)
  expect_s3_class(assess_convergence(m), "gg")
})

test_that("compute_mallows treats observation_frequency properly", {
  dat <- setup_rank_data(potato_visual, observation_frequency = rep(1, nrow(potato_visual)))
  set.seed(123)
  m1 <- compute_mallows(
    data = dat, compute_options = set_compute_options(nmc = 10),
     cl = cl
  )
  dat <- setup_rank_data(potato_visual)
  set.seed(123)
  m2 <- compute_mallows(
    data = dat, compute_options = set_compute_options(nmc = 10),
    cl = cl
  )
  expect_equal(m1, m2)

  # Test with repeated beach preferences
  observation_frequency <- c(2, 1, 4)
  beach_small <- subset(beach_preferences, assessor %in% c(1, 2, 3))

  # Next, we create a new hypothetical beach_preferences dataframe where each
  # assessor is replicated 1-4 times

  beach_pref_rep <- do.call(rbind, lapply(split(beach_small, f = seq_len(nrow(beach_small))), function(dd) {
    ret <- merge(
      dd,
      data.frame(new_assessor = seq_len(observation_frequency[dd$assessor])),
      all = TRUE
    )
    ret$assessor <- paste(ret$assessor, ret$new_assessor, sep = ",")
    ret$new_assessor <- NULL
    ret
  }))

  set.seed(1223)
  dat_small <- setup_rank_data(preferences = beach_small, observation_frequency = observation_frequency)
  dat_rep <- setup_rank_data(preferences = beach_pref_rep)

  model_fit_obs_freq <- compute_mallows(
    data = dat_small,
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    cl = cl
  )

  expect_equal(
    model_fit_obs_freq$rho$value,
    c(7, 1, 14, 9, 15, 10, 2, 8, 5, 11, 4, 6, 13, 3, 12, 7, 1, 14,
      9, 15, 10, 2, 8, 6, 11, 4, 5, 13, 3, 12, 7, 1, 14, 9, 15, 10,
      2, 8, 5, 11, 4, 6, 13, 3, 12, 7, 1, 14, 9, 15, 10, 2, 8, 5, 11,
      4, 6, 13, 3, 12, 7, 1, 14, 9, 15, 10, 2, 8, 5, 11, 4, 6, 13,
      3, 12, 7, 2, 14, 9, 15, 10, 1, 8, 5, 11, 4, 6, 13, 3, 12, 7,
      2, 13, 9, 15, 10, 1, 8, 5, 11, 4, 6, 14, 3, 12, 7, 1, 13, 9,
      15, 10, 2, 8, 5, 11, 4, 6, 14, 3, 12, 7, 1, 13, 9, 15, 11, 2,
      8, 5, 10, 4, 6, 14, 3, 12, 6, 1, 13, 9, 15, 11, 2, 8, 5, 10,
      4, 7, 14, 3, 12, 8, 13, 6, 3, 15, 10, 5, 11, 9, 2, 7, 4, 1, 12,
      14, 8, 13, 6, 3, 15, 11, 5, 10, 9, 2, 7, 4, 1, 12, 14, 8, 13,
      6, 2, 15, 11, 5, 10, 9, 3, 7, 4, 1, 12, 14, 8, 13, 6, 2, 15,
      10, 5, 11, 9, 3, 7, 4, 1, 12, 14, 8, 13, 7, 2, 15, 10, 5, 11,
      9, 3, 6, 4, 1, 12, 14, 8, 13, 7, 2, 14, 10, 5, 11, 9, 3, 6, 4,
      1, 12, 15, 8, 13, 6, 2, 14, 10, 5, 11, 9, 3, 7, 4, 1, 12, 15,
      8, 13, 6, 2, 14, 10, 5, 11, 9, 3, 7, 4, 1, 12, 15, 7, 13, 6,
      2, 14, 10, 5, 11, 9, 3, 8, 4, 1, 12, 15, 7, 13, 5, 2, 14, 10,
      6, 11, 9, 3, 8, 4, 1, 12, 15)
  )

  expect_equal(
    model_fit_obs_freq$alpha$value,
    c(1, 0.993700974746543, 1.0281926554902, 1.0281926554902, 1.09951309221072,
      1.06973509558279, 1.04704509303566, 1.04704509303566, 1.03314038637389,
      1.00089665880574, 1, 1.13001537028802, 1.04557677544508, 1.04579825108138,
      1.04579825108138, 1.07539967807669, 1.07539967807669, 1.10601590595713,
      1.08326163458751, 1.08326163458751)
  )

  # Next for the repeated data.
  set.seed(11)
  model_fit_rep <- compute_mallows(
    data = dat_rep,
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    cl = cl
  )

  expect_equal(
    model_fit_rep$rho$value,
    c(4, 11, 3, 10, 9, 15, 5, 13, 7, 12, 6, 1, 14, 2, 8, 4, 11, 3,
      10, 9, 15, 5, 13, 6, 12, 7, 1, 14, 2, 8, 4, 11, 3, 10, 8, 15,
      5, 13, 6, 12, 7, 1, 14, 2, 9, 4, 11, 3, 10, 8, 15, 5, 14, 6,
      12, 7, 1, 13, 2, 9, 4, 11, 3, 9, 8, 15, 5, 14, 6, 12, 7, 1, 13,
      2, 10, 4, 11, 3, 9, 8, 15, 5, 14, 6, 13, 7, 1, 12, 2, 10, 4,
      11, 3, 9, 8, 15, 5, 14, 6, 13, 7, 2, 12, 1, 10, 4, 11, 3, 9,
      8, 15, 5, 14, 6, 13, 7, 2, 12, 1, 10, 4, 12, 3, 9, 8, 15, 5,
      14, 6, 13, 7, 2, 11, 1, 10, 4, 11, 3, 9, 8, 15, 5, 14, 6, 13,
      7, 2, 12, 1, 10, 2, 3, 10, 13, 15, 12, 8, 1, 9, 11, 6, 4, 7,
      14, 5, 2, 3, 10, 13, 15, 12, 8, 1, 9, 11, 5, 4, 7, 14, 6, 1,
      3, 10, 13, 15, 12, 8, 2, 9, 11, 5, 4, 7, 14, 6, 1, 2, 10, 13,
      15, 12, 8, 3, 9, 11, 5, 4, 7, 14, 6, 1, 2, 9, 13, 15, 12, 8,
      3, 10, 11, 5, 4, 7, 14, 6, 1, 2, 9, 13, 14, 12, 8, 3, 10, 11,
      5, 4, 7, 15, 6, 1, 2, 9, 14, 13, 12, 8, 3, 10, 11, 5, 4, 7, 15,
      6, 1, 2, 9, 15, 13, 12, 8, 3, 10, 11, 5, 4, 7, 14, 6, 1, 2, 9,
      15, 13, 12, 8, 3, 10, 11, 6, 4, 7, 14, 5, 1, 2, 9, 15, 13, 12,
      8, 3, 10, 11, 5, 4, 7, 14, 6)
  )

  expect_equal(
    model_fit_rep$alpha$value,
    c(1, 0.960675800962132, 1.10584778042687, 1.0665805995754, 1.0665805995754,
      1.0665805995754, 1.07184101215446, 0.974563153536915, 1.00698574079831,
      0.869742224461801, 1, 1.0596448436979, 1.0596448436979, 0.878848116239248,
      0.876794712836932, 0.867647585063925, 0.831343269235968, 0.867329621800878,
      0.798324101632821, 0.764940732160401)
  )
})

test_that("compute_mallows() takes initial alpha", {
  mm <- compute_mallows(setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 5),
    initial_values = set_initial_values(alpha_init = c(3, 4)),
    cl = cl
  )
  expect_s3_class(mm, "BayesMallows")
})

parallel::stopCluster(cl)
