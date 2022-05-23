context("SMC new users partial rankings")

# Generate Dataset =======================================

# General ------------------------------------------------
n_items <- dim(sushi_rankings)[2] # Number of items
rho_0 <- seq(from = 1, to = n_items, by = 1) # 'true' consensus ranking
alpha_0 <- 2 # fixed/ 'true' scale parameter
leap_size <- floor(n_items / 5)
metric <- "footrule"

# Generate estimate of Z_n(alpha) ------------------------
alpha_vector <- seq(from = 0, to = 20, by = 1)
iter <- 1e2
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model
# using the estimate partition function
set.seed(994)
logz_estimate <- estimate_partition_function(
  method       = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items      = n_items,
  metric       = metric,
  nmc          = iter,
  degree       = degree
)

# Make this information partially observed over time -----
samples <- sushi_rankings[1:100, ]
samples[samples > 5] <- NA

# Bayesmallows MCMC Results ==============================
nmc <- 2000
bm_mcmc <- compute_mallows(
  samples,
  nmc           = nmc,
  leap_size     = leap_size,
  metric        = metric,
  alpha_prop_sd = 0.15
)
bm_mcmc$burnin <- 1000

# choice items to see in trace plot
items <- sample(1:n_items, 5, replace = FALSE)
items <- sort(items)

test_that("BayesMallows MCMC Results are OK", {
  rho_cp <- compute_consensus(
    model_fit = bm_mcmc,
    type      = "CP",
    burnin    = bm_mcmc$burnin
  )
  rho_map <- compute_consensus(
    model_fit = bm_mcmc,
    type      = "MAP",
    burnin    = bm_mcmc$burnin
  )
  post_rho <- compute_posterior_intervals(bm_mcmc, parameter = "rho")
  post_alpha <- compute_posterior_intervals(bm_mcmc, parameter = "alpha")
  expect_equal(dim(post_rho), c(10, 7))
  expect_equal(dim(rho_cp), c(10, 3))
  expect_equal(dim(rho_map), c(10, 3))
  expect_equal(dim(post_alpha), c(1, 6))
})

# SMC Analysis (alpha unknown) ===========================

mcmc_times <- 5
num_new_obs <- 5
Time <- dim(samples)[1] / num_new_obs
alpha_prop_sd <- 0.5
lambda <- 0.15
alpha_max <- 1e0

test_that("Produces the wrong metric and aug_method error", {
  N <- 5
  expect_error(
    smc_mallows_new_users(
      alpha           = alpha_0,
      type            = "partial_alpha_fixed",
      R_obs           = samples,
      n_items         = n_items,
      metric          = "cayley",
      leap_size       = leap_size,
      N               = N,
      Time            = Time,
      logz_estimate   = logz_estimate,
      mcmc_kernel_app = mcmc_times,
      num_new_obs     = num_new_obs,
      aug_method      = "pseudolikelihood"
    )
  )
  expect_error(
    smc_mallows_new_users(
      R_obs           = samples,
      type            = "partial",
      n_items         = n_items,
      metric          = "cayley",
      leap_size       = leap_size,
      N               = N,
      Time            = Time,
      logz_estimate   = logz_estimate,
      mcmc_kernel_app = mcmc_times,
      num_new_obs     = num_new_obs,
      alpha_prop_sd   = alpha_prop_sd,
      lambda          = lambda,
      alpha_max       = alpha_max,
      aug_method      = "pseudolikelihood"
    )
  )
})
test_that("Runs with unif kernel", {
  N <- 5
  smc_unif_alpha_fixed_unif <- smc_mallows_new_users(
    alpha           = alpha_0,
    type            = "partial_alpha_fixed",
    R_obs           = samples,
    n_items         = n_items,
    metric          = "footrule",
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    aug_method      = "random"
  )
  expect_warning(
    smc_mallows_new_users_partial_alpha_fixed(
      alpha           = alpha_0,
      R_obs           = samples,
      n_items         = n_items,
      metric          = "footrule",
      leap_size       = leap_size,
      N               = N,
      Time            = Time,
      logz_estimate   = logz_estimate,
      mcmc_kernel_app = mcmc_times,
      num_new_obs     = num_new_obs,
      aug_method      = "random"
    ),
    "'smc_mallows_new_users_partial_alpha_fixed' is deprecated."
  )
  expect_is(smc_unif_alpha_fixed_unif, "list")
  expect_equal(length(smc_unif_alpha_fixed_unif), 4)
  expect_equal(dim(smc_unif_alpha_fixed_unif$rho_samples), c(N, 10, 21))
  smc_unif <- smc_mallows_new_users(
    R_obs           = samples,
    type            = "partial",
    n_items         = n_items,
    metric          = "footrule",
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    alpha_prop_sd   = alpha_prop_sd,
    lambda          = lambda,
    alpha_max       = alpha_max,
    aug_method      = "random"
  )
  expect_warning(
    smc_mallows_new_users_partial(
      R_obs           = samples,
      n_items         = n_items,
      metric          = "footrule",
      leap_size       = leap_size,
      N               = N,
      Time            = Time,
      logz_estimate   = logz_estimate,
      mcmc_kernel_app = mcmc_times,
      num_new_obs     = num_new_obs,
      alpha_prop_sd   = alpha_prop_sd,
      lambda          = lambda,
      alpha_max       = alpha_max,
      aug_method      = "random"
    ),
    "'smc_mallows_new_users_partial' is deprecated."
  )
  expect_is(smc_unif, "list")
  expect_equal(length(smc_unif), 4)
  expect_equal(dim(smc_unif$rho_samples), c(N, 10, 21))
  expect_equal(dim(smc_unif$alpha_samples), c(N, 21))

  expect_s3_class(
    plot_alpha_posterior(smc_unif$alpha_samples[, Time + 1], nmc = N, burnin = 2),
    "ggplot"
  )

  expect_s3_class(
    plot_rho_posterior(smc_unif$rho_samples[, , Time + 1], nmc = N, burnin = 2, C = 1),
    "ggplot"
  )
})

test_that("Runs with pseudo kernel", {
  N <- 5
  smc_unif_alpha_fixed_pseudo <- smc_mallows_new_users(
    alpha           = alpha_0,
    type            = "partial_alpha_fixed",
    R_obs           = samples,
    n_items         = n_items,
    metric          = "footrule",
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    aug_method      = "pseudolikelihood"
  )
  expect_is(smc_unif_alpha_fixed_pseudo, "list")
  expect_equal(length(smc_unif_alpha_fixed_pseudo), 4)
  expect_equal(dim(smc_unif_alpha_fixed_pseudo$rho_samples), c(N, 10, 21))
  smc_pseudo <- smc_mallows_new_users(
    R_obs           = samples,
    type            = "partial",
    n_items         = n_items,
    metric          = "footrule",
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    alpha_prop_sd   = alpha_prop_sd,
    lambda          = lambda,
    alpha_max       = alpha_max,
    aug_method      = "pseudolikelihood"
  )
  expect_is(smc_pseudo, "list")
  expect_equal(length(smc_pseudo), 4)
  expect_equal(dim(smc_pseudo$rho_samples), c(N, 10, 21))
  expect_equal(dim(smc_pseudo$alpha_samples), c(N, 21))
})

# Specific example with random aug method ----------------

test_that("Specific example results are OK", {
  N <- 10
  aug_method <- "random"
  set.seed(5482) # necessary for reproducibility of the random aug_method

  test <- smc_mallows_new_users(
    R_obs           = samples,
    type            = "partial",
    n_items         = n_items,
    metric          = metric,
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    alpha_prop_sd   = alpha_prop_sd,
    lambda          = lambda,
    alpha_max       = alpha_max,
    aug_method      = aug_method
  )

  expect_equivalent(
    test$alpha_samples[1:3, 1:4],
    structure(c(0.127095626667142, 1.01987919207187, 0.359836369752884,
                0.252836572866912, 0.707217839795325, 0.337840167678165, 0.123483377876088,
                0.0766661493660047, 0.901749237961603, 0.155759850023329, 0.346042115242859,
                0.0934527121554998), dim = 3:4),
    tolerance = 1e-5
  )

  expect_error(
    compute_rho_consensus(output = test$rho_samples[, , Time + 1], nmc = N, burnin = NULL, C = 1,
                          type = "CP"),
    "Please specify the burnin."
    )

  expect_error(
    compute_rho_consensus(
      output = test$rho_samples[, , Time + 1], nmc = N, burnin = NULL, C = 1,
      type = "MAP"
    ),
    "Please specify the burnin."
  )

  rho_cp <- compute_rho_consensus(
    output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1,
    type = "CP"
  )

  set.seed(545)
  rho_map <- compute_rho_consensus(
    output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1,
    type = "MAP"
  )
  post_rho <- compute_posterior_intervals_rho(
    output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0
  )
  post_alpha <- compute_posterior_intervals_alpha(
    output = test$alpha_samples[, Time + 1], nmc = N, burnin = 0
  )
  expect_equal(dim(rho_cp), c(10, 3))
  expect_equal(dim(rho_map), c(100, 3))
  expect_equal(dim(post_rho), c(10, 7))
  expect_equal(dim(post_alpha), c(1, 6))

  # Test with nonzero burnin
  rho_cp <- compute_rho_consensus(
    output = test$rho_samples[, , Time + 1], nmc = N, burnin = 2, C = 1,
    type = "CP"
  )
  expect_equal(rho_cp$cumprob,
               c(1, 0.625, 0.625, 0.625, 1, 0.5, 0.5, 0.25, 0.125, 0.25))

  rho_map <- compute_rho_consensus(
    output = test$rho_samples[, , Time + 1], nmc = N, burnin = 2, C = 1,
    type = "MAP"
  )
  expect_equal(
    rho_map$probability,
    c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125)
  )

  test_fixed <- smc_mallows_new_users(
    R_obs           = samples,
    type            = "partial_alpha_fixed",
    n_items         = n_items,
    metric          = metric,
    leap_size       = leap_size,
    N               = N,
    Time            = Time,
    logz_estimate   = logz_estimate,
    mcmc_kernel_app = mcmc_times,
    num_new_obs     = num_new_obs,
    aug_method      = aug_method,
    alpha           = alpha_0
  )
  rho_cp_fixed <- compute_rho_consensus(
    output = test_fixed$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1,
    type = "CP"
  )
  set.seed(584)
  rho_map_fixed <- compute_rho_consensus(
    output = test_fixed$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1,
    type = "MAP"
  )
  post_rho_fixed <- compute_posterior_intervals_rho(
    output = test_fixed$rho_samples[, , Time + 1], nmc = N, burnin = 0
  )
  expect_equal(dim(rho_cp_fixed), c(10, 3))
  expect_equal(dim(rho_map_fixed), c(10, 3))
  expect_equal(dim(post_rho_fixed), c(10, 7))

})
