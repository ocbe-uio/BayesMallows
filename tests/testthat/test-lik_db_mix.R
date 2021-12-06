test_that("lik_db_mix works", {
  set.seed(1)
  n_items <- 5
  mydata <- sample_mallows(
    n_samples = 100,
    rho0 = 1:n_items,
    alpha0 = 10,
    metric = "kendall")

  # Compute the likelihood and log-likelihood values under the true model...
  expect_equal(
    sprintf("%.3e", lik_db_mix(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = mydata
  )), "1.434e-74")

  expect_equal(round(lik_db_mix(
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
    sprintf("%.3e", lik_db_mix(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = freq_distr[, 1:n_items],
    obs_freq = freq_distr[, n_items + 1]
  )), "1.434e-74")

  expect_equal(round(lik_db_mix(
    rho = rbind(1:n_items, 1:n_items),
    alpha = c(2 * n_items, 2 * n_items),
    weights = c(0.5, 0.5),
    metric = "kendall",
    rankings = freq_distr[, 1:n_items],
    obs_freq = freq_distr[, n_items + 1],
    log = TRUE
  ), 4), -170.0306)

  expect_error(
    lik_db_mix(
      rho = rbind(1:n_items, 1:n_items),
      alpha = c(2 * n_items, 2 * n_items),
      weights = c(0.5, 0.5),
      metric = "kendall",
      rankings = mydata,
      obs_freq = c(1, 2)
    ))

})
