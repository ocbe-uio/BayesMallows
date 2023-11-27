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
    setup_rank_data(namat, na_action = "augment"),
    "BayesMallowsData"
  )
  expect_s3_class(setup_rank_data(namat), "BayesMallowsData")
  expect_error(
    compute_mallows(
      compute_options = set_compute_options(nmc = 1000, alpha_prop_sd = 1)
    )
  )
  expect_error(
    compute_mallows(
      setup_rank_data(potato_visual),
      compute_options = set_compute_options(nmc = 100, alpha_jump = 102)
    ),
    "nmc must be strictly larger than alpha_jump"
  )
  expect_error(
    compute_mallows(
      data = setup_rank_data(potato_visual),
      priors = set_priors(lambda = 0)
    ),
    "lambda must be a strictly positive number of length one"
  )
  expect_error(
    compute_mallows(
      data = setup_rank_data(potato_visual),
      priors = set_priors(lambda = -10.2)
    ),
    "lambda must be a strictly positive number of length one"
  )
  expect_error(
    compute_mallows(
      data = setup_rank_data(potato_visual),
      compute_options = set_compute_options(nmc = 100, rho_thinning = 200)
    ),
    "nmc must be strictly larger than rho_thinning"
  )
  expect_error(
    compute_mallows(
      data = setup_rank_data(potato_visual),
      compute_options = set_compute_options(nmc = 100, aug_thinning = 200)
    ),
    "nmc must be strictly larger than aug_thinning"
  )
  expect_error(
    compute_mallows(
      data = setup_rank_data(potato_visual),
      compute_options = set_compute_options(nmc = -100)
    ),
    "nmc must be a positive integer"
  )
})

test_that("rho_init is properly validated", {
  m <- potato_visual
  expect_error(
    compute_mallows(
      data = setup_rank_data(m),
      initial_values = set_initial_values(rho_init = 1:(ncol(m) - 1))
    ),
    "initial value rho must have one value per item"
  )

  expect_error(
    compute_mallows(
      data = setup_rank_data(m),
      initial_values = set_initial_values(rho_init = c(potato_true_ranking[-1], 22))
    ),
    "rho_init must be a proper permutation"
  )

  expect_error(
    compute_mallows(
      data = setup_rank_data(m),
      initial_values = set_initial_values(rho_init = c(NA_real_, 2:ncol(m)))
    ),
    "rho_init cannot have missing values"
  )

  expect_error(
    setup_rank_data(m, observation_frequency = -1),
    "observation_frequency must be a vector of strictly positive numbers"
  )

  expect_error(
    setup_rank_data(m, observation_frequency = 1),
    "observation_frequency must be of same length as the number of rows in rankings"
  )

  expect_error(
    setup_rank_data(m, observation_frequency = 1:11),
    "observation_frequency must be of same length as the number of rows in rankings"
  )
})

test_that("inconsistent rankings are discovered", {
  expect_error(
    setup_rank_data(matrix(c(
      1, 2, -3,
      1, 2, 3
    ), nrow = 2, byrow = TRUE)),
    "invalid permutations provided in rankings matrix"
  )

  expect_error(
    setup_rank_data(matrix(c(
      1, 2, 2,
      1, 2, 3
    ), nrow = 2, byrow = TRUE)),
    "invalid permutations provided in rankings matrix"
  )
})


test_that("compute_mallows error model works", {
  preferences <- data.frame(
    assessor = c(1, 1, 2, 2),
    bottom_item = c(1, 2, 1, 2),
    top_item = c(2, 1, 2, 3)
  )
  expect_error(
    compute_mallows(
      setup_rank_data(preferences = preferences),
      compute_options = set_compute_options(nmc = 10)
    ),
    "Intransitive pairwise comparisons. Please specify an error model."
  )
  expect_message(
    compute_mallows(
      data = setup_rank_data(preferences = preferences),
      model_options = set_model_options(error_model = "bernoulli"),
      compute_options = set_compute_options(nmc = 10)
    ),
    "Preferences are intransitive."
  )
})


test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(setup_rank_data(mat),
    compute_options = set_compute_options(nmc = 30)
  )
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")
})


test_that("compute_mallows runs with the right distances", {
  dat <- setup_rank_data(potato_visual)
  for (metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")) {
    expect_s3_class(
      compute_mallows(dat,
        model_options = set_model_options(metric = metric),
        compute_options = set_compute_options(nmc = 3)
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

  expect_s3_class(
    compute_mallows(
      setup_rank_data(preferences = m),
      compute_options = set_compute_options(nmc = 20)
    ), "BayesMallows"
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

  m <- compute_mallows(setup_rank_data(R_partial2))
  expect_s3_class(assess_convergence(m), "gg")
})

test_that("compute_mallows treats observation_frequency properly", {
  set.seed(2233)
  m1 <- compute_mallows(
    setup_rank_data(potato_visual,
      observation_frequency = rep(1, nrow(potato_visual))
    )
  )
  set.seed(2233)
  m2 <- compute_mallows(
    setup_rank_data(potato_visual)
  )
  expect_equal(m1, m2)

  # Test with repeated beach preferences
  observation_frequency <- c(2, 1, 4)

  beach_small <- subset(beach_preferences, assessor %in% c(1, 2, 3))

  # Next, we create a new hypthetical beach_preferences dataframe where each
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


})

