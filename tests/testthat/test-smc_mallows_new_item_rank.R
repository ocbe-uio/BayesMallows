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
    ),
    "Combined choice of metric and aug_method is incompatible"
  )
  expect_error(
    smc_mallows_new_item_rank(
      n_items = n_items, R_obs = sample_dataset,
      metric = "cayley", leap_size = leap_size, N = N, Time = Time,
      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
      alpha_prop_sd = alpha_prop_sd, lambda = lambda,
      alpha_max = alpha_max, aug_method = "pseudolikelihood"
    ),
    "Combined choice of metric and aug_method is incompatible"
  )
})

test_that("smc_mallows_new_item_rank_alpha_fixed is deprecated", {
  expect_warning(
    smc_mallows_new_item_rank_alpha_fixed(
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
  expect_is(smc_unif_alpha_fixed_unif, "SMCMallows")
  expect_length(smc_unif_alpha_fixed_unif, 3)
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
  expect_is(smc_unif, "SMCMallows")
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
  expect_is(smc_unif_alpha_fixed_unif, "SMCMallows")
  expect_length(smc_unif_alpha_fixed_unif, 3)
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
  expect_is(smc_unif, "SMCMallows")
  expect_length(smc_unif, 4)
  expect_equal(dim(smc_unif$rho_samples), c(N, 6, 31))
  expect_equal(dim(smc_unif$alpha_samples), c(N, 31))
  expect_equal(smc_unif$augmented_rankings[, , 10],
               structure(c(1, 2, 3, 1, 3, 3, 2, 1, 1, 1, 2, 4, 2, 5, 1, 1, 3,
                           3, 2, 2, 3, 5, 4, 3, 5, 2, 5, 2, 6, 3, 4, 6, 1, 4, 6, 4, 4, 4,
                           5, 4, 6, 1, 5, 2, 4, 5, 1, 5, 3, 5, 5, 3, 6, 6, 2, 6, 6, 6, 4,
                           6), dim = c(10L, 6L)))
  expect_equal(smc_unif$augmented_rankings[, , 17],
                    structure(c(1, 2, 3, 1, 3, 3, 2, 1, 1, 1, 2, 4, 2, 5, 1, 1, 3,
                                3, 2, 2, 3, 5, 4, 3, 5, 2, 5, 2, 6, 3, 4, 6, 1, 4, 6, 4, 4, 4,
                                5, 4, 6, 1, 5, 2, 4, 5, 1, 5, 3, 5, 5, 3, 6, 6, 2, 6, 6, 6, 4,
                                6), dim = c(10L, 6L)))
  expect_equal(
    smc_unif$rho_samples[, , 7],
    structure(
      c(
        4, 6, 6, 6, 3, 6, 6, 6, 6, 6, 4, 6, 6, 5, 6, 5, 5, 5, 5, 5, 5, 3, 5, 5,
        4, 5, 3, 3, 4, 5, 3, 5, 5, 6, 4, 3, 4, 4, 3, 4, 3, 1, 1, 1, 1, 1, 1, 1,
        1, 2, 1, 4, 1, 1, 1, 2, 2, 2, 4, 2, 6, 4, 4, 4, 5, 4, 4, 4, 5, 3, 6, 2,
        4, 4, 5, 4, 6, 6, 6, 6, 2, 5, 2, 2, 6, 2, 5, 5, 3, 1, 5, 3, 2, 2, 3, 6,
        3, 3, 1, 3, 1, 2, 3, 3, 2, 3, 2, 2, 2, 4, 2, 1, 3, 3, 2, 1, 1, 1, 2, 1
      ),
      dim = c(20L, 6L)
    )
  )
  expect_equal(
    smc_unif$alpha_samples[, 10],
    c(
      0.1455621, 0.1945882, 0.0784438, 0.1732896, 0.09138699, 0.2689943,
      0.09694153, 0.1892184, 0.505403, 0.07356979, 0.7224048, 0.5333023,
      0.1870632, 0.3981105, 0.1357151, 0.318864, 0.1283664, 0.3174649,
      0.6529184, 0.08086917
    ),
    tol = 1e-6
  )
})
