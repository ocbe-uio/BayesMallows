test_that("get_mallows_loglik works", {
  set.seed(1)
  n_items <- 5
  mydata <- sample_mallows(
    n_samples = 100,
    rho0 = 1:n_items,
    alpha0 = 10,
    metric = "kendall"
  )

  # Compute the likelihood and log-likelihood values under the true model...
  # Two mixture components
  expect_equal(
    lapply(
      c("ulam", "footrule", "spearman", "kendall", "cayley", "hamming"),
      function(m){
        get_mallows_loglik(
          rho = rbind(1:n_items, 1:n_items),
          alpha = c(2 * n_items, 2 * n_items),
          weights = c(0.5, 0.5),
          metric = m,
          rankings = mydata,
          log = TRUE
        )
      }),
    list(-241.809448010519, -231.458744596071, -283.167543856933,
         -170.030563541214, -223.987863416978, -224.300341956488))

  # A single component
  # Gives exactly the same values as the previous, due to the weights
  expect_equal(
    lapply(
      c("ulam", "footrule", "spearman", "kendall", "cayley", "hamming"),
      function(m){
        get_mallows_loglik(
          rho = 1:n_items,
          alpha = 2 * n_items,
          weights = 1,
          metric = m,
          rankings = mydata,
          log = TRUE
        )
      }),
    list(-241.809448010519, -231.458744596071, -283.167543856933,
         -170.030563541214, -223.987863416978, -224.300341956488))

  expect_equal(
    sprintf("%.3e", get_mallows_loglik(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = mydata,
      log = FALSE
    )), "1.434e-74"
  )

  expect_equal(round(get_mallows_loglik(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = mydata,
    log = TRUE
  ), 4), -170.0306)

  # or equivalently, by using the frequency distribution
  freq_distr <- rank_freq_distr(mydata)
  expect_equal(
    sprintf("%.3e", get_mallows_loglik(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = freq_distr[, 1:n_items],
      obs_freq = freq_distr[, n_items + 1],
      log = FALSE
    )), "1.434e-74"
  )

  expect_equal(round(get_mallows_loglik(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = freq_distr[, 1:n_items],
    obs_freq = freq_distr[, n_items + 1],
    log = TRUE
  ), 4), -170.0306)

  expect_error(
    get_mallows_loglik(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = mydata,
      obs_freq = c(1, 2)
    ),
    "obs_freq must be of same length as the number of rows in rankings"
  )

  expect_error(
    get_mallows_loglik(
      rho = 1:100,
      alpha = 1,
      weights = 1,
      metric = "spearman",
      rankings = do.call(rbind, replicate(3, list(1:100), simplify = "list"))
    ),
    "Given number of items currently not available for the specified metric"
  )

  expect_warning(
    lik_db_mix(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = freq_distr[, 1:n_items],
      obs_freq = freq_distr[, n_items + 1],
      log = TRUE
    ),
    "lik_db_mix is deprecated, and the log argument now defaults to TRUE."
  )
})
