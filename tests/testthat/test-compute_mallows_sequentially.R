test_that("compute_mallows_sequentially works", {
  set.seed(345)
  data <- lapply(seq_len(nrow(potato_visual)), function(i) {
    setup_rank_data(potato_visual[i, ])
  })

  initial_values <- sample_prior(
    n = 200, n_items = 20,
    priors = set_priors(gamma = 3, lambda = .1)
  )

  mod <- compute_mallows_sequentially(
    data = data,
    initial_values = initial_values,
    smc_options = set_smc_options(n_particles = 200, mcmc_steps = 20)
  )

  expect_equal(
    apply(mod$alpha_samples, 2, mean),
    c(
      2.71999372795106, 2.19488562724823, 2.31651832681418, 3.38677205663342,
      3.94905442011209, 6.37972319012761, 8.38459498953842, 9.37799520951236,
      9.88211129583396, 10.8952555605406, 10.7448627204384, 10.9013851034043
    )
  )

  expect_equal(
    mod$rho_samples[4, c(3, 9), c(1, 9, 10)],
    matrix(c(9, 18, 11, 8, 16, 16), ncol = 3)
  )

  expect_equal(get_acceptance_ratios(mod)$alpha_acceptance[[5]], .80025)
  expect_equal(get_acceptance_ratios(mod)$rho_acceptance[[7]], .3495)

  expect_error(
    compute_mallows_sequentially(potato_visual, initial_values),
    "data must be a list of BayesMallowsData objects."
  )

  expect_error(
    compute_mallows_sequentially(setup_rank_data(potato_visual), initial_values),
    "data must be a list of BayesMallowsData objects."
  )
})

test_that("compute_mallows_sequentially works with partial rankings", {
  set.seed(345)
  dat <- potato_visual
  dat[dat > 15] <- NA

  data <- lapply(seq_len(nrow(dat)), function(i) {
    setup_rank_data(dat[i, ])
  })

  initial_values <- sample_prior(
    n = 200, n_items = 20,
    priors = set_priors(gamma = 3, lambda = .1)
  )

  mod <- compute_mallows_sequentially(
    data = data,
    initial_values = initial_values,
    smc_options = set_smc_options(n_particles = 200, mcmc_steps = 20)
  )

  expect_equal(
    apply(mod$alpha_samples, 2, mean),
    c(
      2.87674720640401, 2.22027499892846, 2.53499158193073, 3.28424790116212,
      3.86836385208775, 4.57947369845994, 5.75771795091753, 7.5070197872037,
      8.63830524611693, 9.73845592788363, 10.2889617544371, 10.4747448422823
    )
  )
})
