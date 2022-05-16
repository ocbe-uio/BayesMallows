context("SMC new user and item rank combined")

# a simpler example to test ====================================================
set.seed(101)
Time <- dim(sample_dataset)[3]

# General ======================================================================
n_items <- dim(sample_dataset)[2] # Number of items
rho_0 <- seq(from = 1, to = n_items, by = 1) # 'true' consensus ranking
alpha_0 <- 2 # fixed/ 'true' scale parameter
leap_size <- floor(n_items / 5)
metric <- "footrule"

# Generate estimate of Z_n(alpha) ==============================================
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e2
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model
# using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)

mcmc_kernel_app <- 5
N <- 20
alpha_prop_sd <- 0.5
lambda <- 0.15
alpha_max <- 1e6

test_that("Produces the wrong metric and aug_method error", {
  expect_error(
    smc_mallows_new_item_rank(
      alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
      metric = "cayley", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "pseudolikelihood",
      alpha_fixed = TRUE
    )
  )
  expect_error(
    smc_mallows_new_item_rank(
      n_items = n_items, R_obs = sample_dataset,
      metric = "cayley", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "pseudolikelihood"
    )
  )
})

test_that("smc_mallows_new_item_rank_alpha_fixed is deprecated", {
  expect_warning(smc_mallows_new_item_rank_alpha_fixed(
    alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
    metric = "footrule", leap_size = leap_size, N = N, Time = Time,
    logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
    alpha_prop_sd = alpha_prop_sd, lambda = lambda,
    alpha_max = alpha_max, aug_method = "random"
  ),
  "'smc_mallows_new_item_rank_alpha_fixed' is deprecated."
  )
})

test_that("Runs with unif kernel", {
  set.seed(0)
  smc_unif_alpha_fixed_unif <- suppressMessages(
    smc_mallows_new_item_rank(
      alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
      metric = "footrule", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "random", alpha_fixed = TRUE
    )
  )
  expect_is(smc_unif_alpha_fixed_unif, "list")
  expect_length(smc_unif_alpha_fixed_unif, 4)
  expect_equal(dim(smc_unif_alpha_fixed_unif$rho_samples), c(N, 6, 31))
  expect_equal(
    smc_unif_alpha_fixed_unif$augmented_rankings[c(4, 7), 5, c(13, 19)],
    structure(c(2, 1, 2, 1), dim = c(2L, 2L))
  )
  set.seed(2)
  smc_unif <- suppressMessages(
    smc_mallows_new_item_rank(
      n_items = n_items, R_obs = sample_dataset,
      metric = "footrule", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "random"
    )
  )
  expect_is(smc_unif, "list")
  expect_length(smc_unif, 4)
  expect_equal(dim(smc_unif$rho_samples), c(N, 6, 31))
  expect_equal(dim(smc_unif$alpha_samples), c(N, 31))
})

test_that("Runs with pseudo kernel", {
  smc_unif_alpha_fixed_unif <- suppressMessages(
    smc_mallows_new_item_rank(
      alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
      metric = "footrule", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "pseudolikelihood",
      alpha_fixed = TRUE
    )
  )
  expect_is(smc_unif_alpha_fixed_unif, "list")
  expect_length(smc_unif_alpha_fixed_unif, 4)
  expect_equal(dim(smc_unif_alpha_fixed_unif$rho_samples), c(N, 6, 31))

  set.seed(1)
  smc_unif <- suppressMessages(
    smc_mallows_new_item_rank(
      n_items = n_items, R_obs = sample_dataset,
      metric = "footrule", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "pseudolikelihood"
    )
  )
  expect_is(smc_unif, "list")
  expect_length(smc_unif, 4)
  expect_equal(dim(smc_unif$rho_samples), c(N, 6, 31))
  expect_equal(dim(smc_unif$alpha_samples), c(N, 31))
  expect_equal(smc_unif$augmented_rankings[,,10],
               structure(c(1, 2, 3, 1, 3, 3, 2, 1, 1, 1, 2, 4, 2, 5, 1, 1, 3,
                           3, 2, 2, 3, 5, 4, 3, 5, 2, 5, 2, 6, 3, 4, 6, 1, 4, 6, 4, 4, 4,
                           5, 4, 6, 1, 5, 2, 4, 5, 1, 5, 3, 5, 5, 3, 6, 6, 2, 6, 6, 6, 4,
                           6), dim = c(10L, 6L)))
  expect_equal(smc_unif$augmented_rankings[,,17],
                    structure(c(1, 2, 3, 1, 3, 3, 2, 1, 1, 1, 2, 4, 2, 5, 1, 1, 3,
                                3, 2, 2, 3, 5, 4, 3, 5, 2, 5, 2, 6, 3, 4, 6, 1, 4, 6, 4, 4, 4,
                                5, 4, 6, 1, 5, 2, 4, 5, 1, 5, 3, 5, 5, 3, 6, 6, 2, 6, 6, 6, 4,
                                6), dim = c(10L, 6L)))
  expect_equal(
    smc_unif$rho_samples[,,7],
    structure(c(2, 1, 1, 3, 3, 1, 1, 2, 1, 1, 1, 3, 2, 3, 1, 2, 1,
                3, 2, 1, 1, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2,
                1, 2, 4, 4, 4, 5, 5, 4, 3, 6, 3, 5, 5, 4, 4, 4, 4, 6, 5, 4, 6,
                3, 5, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 5, 6, 5, 5, 5, 6, 6, 5, 6,
                3, 3, 3, 1, 1, 3, 4, 1, 5, 4, 4, 1, 5, 6, 3, 4, 4, 1, 4, 4, 6,
                5, 5, 4, 4, 5, 5, 4, 4, 3, 3, 6, 3, 2, 6, 3, 3, 5, 3, 5), dim = c(20L,
                                                                                  6L))
  )
  expect_equal(
    smc_unif$alpha_samples[,10],
    c(0.602219711983957, 0.864155419159704, 0.19987176844615, 0.212150389835116,
      0.215010632154594, 0.214553410823427, 2.44508260024763, 0.588951840125527,
      0.432062764370222, 0.698030157143868, 0.284210349887872, 0.202150647732622,
      0.569533090734916, 0.1463613773319, 0.253695697680041, 0.863234080153403,
      0.545039808362408, 0.201199396612719, 2.10350905472693, 0.371588907112214
    )
  )
})
