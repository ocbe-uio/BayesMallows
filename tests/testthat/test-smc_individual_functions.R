context("SMC individual functions")

rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6

test_that("get_exponent_sum() works as expected", {
  set.seed(101)
  loglik <- get_exponent_sum(
    alpha = alpha, rho = t(rho), n_items = length(rho), rankings = t(rho),
    metric = metric
  )
  expect_equal(loglik, 0)

  rankings <- sample_mallows(
    rho0 = rho, alpha0 = alpha, n_samples = 10,
    burnin = 1000, thinning = 500
  )
  loglik <- get_exponent_sum(
    alpha = alpha, rho = rho, n_items = n_items, rankings = rankings,
    metric = metric
  )
  expect_equivalent(loglik, -22.6667, tol = 1e-4)
})

test_that("smc_metropolis_hastings_rho() works as expected", {
  set.seed(101)
  # This functions uses get_mallows_log_lik and leap_and_shift_probs
  # so if the checks match in those worker functions then it is very likely
  # that this function will return the correct outputs.
  rankings <- sample_mallows(
    rho0 = rho, alpha0 = alpha, n_samples = 10,
    burnin = 1000, thinning = 500
  )

  # you can confirm the print statements inside the metropolis_hastings_rho
  # match get_exponent_sum and leap_and_shift_probs
  test_1 <- metropolis_hastings_rho(
    alpha = alpha, n_items = n_items, rankings = t(rho), metric = metric,
    rho = rho, leap_size = 1
  )
  dist_1 <- get_rank_distance(rho, test_1, metric = "ulam")
  expect_equal(test_1, as.matrix(c(1, 2, 3, 5, 4, 6)))
  # if rho != rho_prime, then it should have a ulam distance of 1
  # if rho == rho_prime, then it should have ulam distance of 0
  expect_equal(dist_1, 1)

  test_2 <- metropolis_hastings_rho(
    alpha = alpha, n_items = n_items, rankings = t(rho), metric = metric,
    rho = rho, leap_size = 2
  )
  dist_2 <- get_rank_distance(rho, test_2, metric = "ulam")
  expect_equal(test_2, as.matrix(c(1, 2, 3, 4, 5, 6)))
  expect_equal(dist_2, 0)

  test_3 <- metropolis_hastings_rho(
    alpha = alpha, n_items = n_items, rankings = t(rho), metric = metric,
    rho = rho, leap_size = 3
  )
  dist_3 <- get_rank_distance(rho, test_3, metric = "ulam")
  expect_equal(test_3, as.matrix(c(1, 2, 3, 4, 5, 6)))
  expect_equal(dist_3, 0)

  # we have a ranking data set containing 10 rankings over 6 items
  test_4 <- metropolis_hastings_rho(
    alpha = alpha, n_items = n_items, rankings = rankings, metric = metric,
    rho = rho, leap_size = 1
  )
  dist_4 <- get_rank_distance(rho, test_4, metric = "ulam")
  expect_equal(test_4, as.matrix(c(1, 2, 3, 4, 5, 6)))
  expect_equal(dist_4, 0)
})

test_that("smc_leap_and_shift_probs() works as expected", {
  set.seed(101)
  n_items <- length(rho)

  # leap_size has a possible range, the BayesMallows papers suggest
  # leap_size = floor(n_items/5) but the leap_size can be up to n_items/2.
  # Note that leap_size must be integered valued.

  # if leap_size = 1, then forwards_prob = backwards_prob
  test_1 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 1)
  expect_equal(test_1$rho_prime, as.matrix(c(1, 3, 2, 4, 5, 6)))
  expect_equivalent(test_1$forwards_prob, 0.1666667, tol = 1e-6)
  expect_equivalent(test_1$backwards_prob, 0.1666667, tol = 1e-6)

  # if rho != rho_prime, then it should have a ulam distance of 1
  # if rho == rho_prime, then it should have ulam distance of 0
  dist_1 <- get_rank_distance(rho, test_1$rho_prime, metric = "ulam")
  expect_equal(dist_1, 1)

  test_2 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 2)
  expect_equal(test_2$rho_prime, as.matrix(c(1, 2, 3, 5, 4, 6)))
  expect_equivalent(test_2$forwards_prob, 0.09722222, tol = 1e-4)
  expect_equivalent(test_2$backwards_prob, 0.09722222, tol = 1e-4)

  dist_2 <- get_rank_distance(rho, test_2$rho_prime, metric = "ulam")
  expect_equal(dist_2, 1)

  test_3 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 3)
  expect_equal(test_3$rho_prime, as.matrix(c(1, 3, 2, 4, 5, 6)))
  expect_equivalent(test_3$forwards_prob, 0.075, tol = 1e-3)
  expect_equivalent(test_3$backwards_prob, 0.075, tol = 1e-3)

  dist_3 <- get_rank_distance(rho, test_3$rho_prime, metric = "ulam")
  expect_equal(dist_3, 1)
})

# ======================================================== #
# unit test for metropolis_hastings_alpha.R                #
# ======================================================== #

metropolis_hastings_alpha_old <- function(alpha, n_items, rankings, metric, rho, logz_estimate) {
  exp_alpha_prime <- rlnorm(1, mean = alpha, sd = 0.15) # 1
  alpha_prime <- log(exp_alpha_prime)

  # evaluate the log-likelihood with current rankings
  mallows_loglik_prop <- get_exponent_sum(
    alpha = (alpha_prime - alpha), rho = rho, n = n_items,
    rankings = rankings, metric = metric
  )

  # evaluate the log estimate of the partition function
  # for a particular value of alpha
  logz_alpha <- get_partition_function(
    n_items = n_items, alpha = alpha, logz_estimate = logz_estimate,
    metric = metric
  )
  logz_alpha_prime <- get_partition_function(
    n_items = n_items, alpha = alpha_prime, logz_estimate = logz_estimate,
    metric = metric
  )

  n_users <- length(rankings) / n_items

  loga <- n_users * (logz_alpha - logz_alpha_prime) +
    dexp(alpha_prime, log = TRUE) - dexp(alpha, log = TRUE) +
    alpha_prime - alpha + mallows_loglik_prop

  # determine whether to accept or reject proposed rho and
  # return now consensus ranking
  p <- runif(1, min = 0, max = 1)
  if (log(p) <= loga) {
    return(alpha_prime)
  } else {
    return(alpha)
  }
}

set.seed(101)
rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6
rankings <- sample_mallows(
  rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
)
alpha_vector <- seq(from = 0, to = 20, by = 1)
iter <- 1e4
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model
# using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling", alpha_vector = alpha_vector,
  n_items = n_items, metric = "footrule", nmc = iter, degree = degree
)
set.seed(101)
test_1_a <- metropolis_hastings_alpha_old(alpha, n_items, rankings, metric, rho, logz_estimate)
test_1_b <- metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate, alpha_prop_sd = 0.5,
  lambda = 0.1, alpha_max = 20, metric
)
set.seed(101)
test_2_a <- metropolis_hastings_alpha_old(
  alpha, n_items, rankings, metric, rho, logz_estimate
)
test_2_b <- metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate, alpha_prop_sd = 0.15,
  lambda = 0.1, alpha_max = 20, metric
)
set.seed(101)
test_3_b <- metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate, alpha_prop_sd = 0.5,
  lambda = 0.15, alpha_max = 20, metric
)
set.seed(101)
test_4_b <- metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate, alpha_prop_sd = 0.15,
  lambda = 0.15, alpha_max = 20, metric
)

test_that("metropolis_hastings_alpha() works as expected", {
  expect_equivalent(test_1_a, 1.951095, tol = 1e-5)
  expect_equivalent(test_1_b, 2.450351, tol = 1e-5)
  expect_equivalent(test_2_a, 1.951095, tol = 1e-5)
  expect_equivalent(test_2_b, 2.125639, tol = 1e-5)
  expect_equivalent(test_3_b, 2)
  expect_equivalent(test_4_b, 1.904542, tol = 1e-5)
})


test_that("leap_and_shift_probs does not propose current ranking", {
  set.seed(12)
  count <- Reduce(`+`, lapply(list(1:4, 1:10), function(rho) {
    n_items <- length(rho)
    count <- 0
    for (i in 1:20) {
      val <- leap_and_shift_probs(rho, 1, n_items)
      if (all(val$rho_prime == rho)) {
        count <- count + 1
      }
      val <- leap_and_shift_probs(rho, 2, n_items)
      if (all(val$rho_prime == rho)) {
        count <- count + 1
      }
    }
    count
  }))
  expect_equal(count, 0)

})
