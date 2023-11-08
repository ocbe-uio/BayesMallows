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
      model = set_model_options(error_model = "bernoulli"),
      compute_options = set_compute_options(nmc = 10)
    ),
    "Preferences are intransitive."
  )
})

test_that("compute_mallows with single missing value works", {
  dd <- potato_visual
  dd[1, 1] <- NA
  dd[2, 3] <- NA
  m <- compute_mallows(
    setup_rank_data(dd),
    compute_options = set_compute_options(nmc = 4), seed = 123L
  )
  expect_equal(
    m$alpha,
    structure(list(
      chain = structure(c(1L, 1L, 1L, 1L), levels = "1", class = "factor"),
      cluster = structure(c(1L, 1L, 1L, 1L), levels = "Cluster 1", class = "factor"),
      iteration = c(1, 2, 3, 4), value = c(
        1, 0.986228529947352,
        0.834184330130122, 0.81366346172066
      )
    ), row.names = c(
      NA,
      -4L
    ), class = "data.frame")
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
        model = set_model_options(metric = metric),
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
  m1 <- compute_mallows(
    setup_rank_data(potato_visual,
      observation_frequency = rep(1, nrow(potato_visual))
    ),
    seed = 2233
  )
  m2 <- compute_mallows(
    setup_rank_data(potato_visual),
    seed = 2233
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

  # We generate the initial rankings for the repeated and the "unrepeated"
  # data
  set.seed(1223)
  dat <- setup_rank_data(preferences = beach_small, observation_frequency = observation_frequency)
  dat_rep <- setup_rank_data(preferences = beach_pref_rep)

  model_fit_obs_freq <- compute_mallows(
    data = dat,
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    seed = 3344L
  )

  expect_equal(
    model_fit_obs_freq$rho$value,
    c(
      15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15, 4, 12,
      13, 8, 6, 2, 14, 3, 10, 11, 1, 9, 7, 5, 15, 4, 12, 13, 7, 6,
      2, 14, 3, 10, 11, 1, 9, 8, 5, 15, 4, 12, 13, 8, 6, 2, 14, 3,
      10, 11, 1, 9, 7, 5, 15, 4, 12, 13, 8, 6, 2, 14, 3, 10, 11, 1,
      9, 7, 5, 15, 5, 12, 13, 8, 6, 2, 14, 3, 10, 11, 1, 9, 7, 4, 15,
      5, 12, 13, 8, 6, 2, 14, 3, 9, 11, 1, 10, 7, 4, 15, 5, 12, 13,
      8, 6, 2, 14, 4, 9, 11, 1, 10, 7, 3, 14, 5, 12, 13, 8, 6, 2, 15,
      4, 9, 11, 1, 10, 7, 3, 14, 5, 12, 13, 8, 6, 2, 15, 3, 9, 11,
      1, 10, 7, 4
    )
  )

  expect_equal(
    model_fit_obs_freq$alpha$value,
    c(
      1, 1.07422639598934, 0.867111587670556, 0.870763091397458,
      0.814708756686028, 0.763045119115647, 0.775326303198559, 0.695442598519739,
      0.67554927571174, 0.600167230862289
    )
  )

  # Next for the repeated data.
  model_fit_rep <- compute_mallows(
    data = dat_rep,
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE),
    seed = 3344L
  )


  expect_equal(
    model_fit_rep$rho$value,
    c(
      15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15, 4, 12,
      13, 8, 6, 3, 14, 2, 11, 10, 1, 9, 7, 5, 15, 4, 11, 13, 8, 6,
      3, 14, 2, 12, 10, 1, 9, 7, 5, 15, 4, 11, 13, 8, 5, 3, 14, 2,
      12, 10, 1, 9, 7, 6, 15, 4, 11, 13, 8, 5, 3, 14, 1, 12, 10, 2,
      9, 7, 6, 14, 4, 11, 13, 8, 5, 3, 15, 1, 12, 10, 2, 9, 7, 6, 14,
      4, 11, 13, 8, 5, 3, 15, 1, 12, 9, 2, 10, 7, 6, 14, 4, 10, 13,
      8, 5, 3, 15, 1, 12, 9, 2, 11, 7, 6, 14, 4, 10, 13, 9, 5, 3, 15,
      1, 12, 8, 2, 11, 7, 6, 14, 4, 10, 13, 9, 5, 3, 15, 1, 12, 8,
      2, 11, 7, 6
    )
  )

  expect_equal(
    model_fit_rep$alpha$value,
    c(
      1, 1.09844540139283, 0.991042262448225, 0.902332143366064,
      1.00941552183767, 1.05591659293716, 1.2025606578257, 1.16106661791628,
      1.25137708104917, 1.0752293035996
    )
  )
})
