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
      function(m) {
        get_mallows_loglik(
          rho = rbind(1:n_items, 1:n_items),
          alpha = c(2 * n_items, 2 * n_items),
          weights = c(0.5, 0.5),
          metric = m,
          rankings = mydata,
          log = TRUE
        )
      }
    ),
    list(
      -241.809448010519, -231.458744596071, -283.167543856933,
      -170.030563541214, -223.987863416978, -224.300341956488
    )
  )

  # A single component
  # Gives exactly the same values as the previous, due to the weights
  expect_equal(
    lapply(
      c("ulam", "footrule", "spearman", "kendall", "cayley", "hamming"),
      function(m) {
        get_mallows_loglik(
          rho = 1:n_items,
          alpha = 2 * n_items,
          weights = 1,
          metric = m,
          rankings = mydata,
          log = TRUE
        )
      }
    ),
    list(
      -241.809448010519, -231.458744596071, -283.167543856933,
      -170.030563541214, -223.987863416978, -224.300341956488
    )
  )

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
  freq_distr <- compute_observation_frequency(mydata)
  expect_equal(
    sprintf("%.3e", get_mallows_loglik(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = freq_distr[, 1:n_items],
      observation_frequency = freq_distr[, n_items + 1],
      log = FALSE
    )), "1.434e-74"
  )

  expect_equal(round(get_mallows_loglik(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = freq_distr[, 1:n_items],
    observation_frequency = freq_distr[, n_items + 1],
    log = TRUE
  ), 4), -170.0306)

  expect_error(
    get_mallows_loglik(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = mydata,
      observation_frequency = c(1, 2)
    ),
    "observation_frequency must be of same length as the number of rows in rankings"
  )

  expect_error(
    get_mallows_loglik(
      rho = 1:100,
      alpha = 1,
      weights = 1,
      metric = "spearman",
      rankings = do.call(rbind, replicate(3, list(1:100), simplify = "list"))
    ),
    "Partition function not available."
  )
})

test_that("get_mallows_loglik is correct with clusters", {
  rankings <- R <- potato_visual
  n_items <- ncol(R)
  N <- nrow(R)

  rho <- rbind(potato_true_ranking,1:20)
  alpha <- c(2.5,1)
  weights <- c(0.2,0.8)

  expect_equal(
    get_mallows_loglik(rho[1,],alpha[1],1,rankings = R,metric='spearman'),
    -279.763590378285)

  expect_equal(
    get_mallows_loglik(rho,alpha,weights,rankings = R,metric='spearman'),
    -299.076845327494)
})
