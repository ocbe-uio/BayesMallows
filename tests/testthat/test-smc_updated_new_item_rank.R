# set-up

# new item rank for each user (fewer things)
example_dataset <- sushi_rankings[1:100, ]
n_users <- dim(example_dataset)[1]
n_items <- dim(example_dataset)[2]
test_dataset <- array(0, c(n_users, n_items, (n_items / 2 + 1)))
test_dataset[, , (n_items / 2 + 1)] <- example_dataset
tt <- 0
for (ii in (n_items - 1):(n_items / 2)) {
  tt <- tt + 1

  # set n_users line with one more NA
  example_dataset[example_dataset > ii] <- NA

  # set as new time stamp
  test_dataset[, , ((n_items / 2 + 1) - tt)] <- example_dataset
}


metric <- "footrule"
leap_size <- floor(n_items / 5)

cardinalities <- prepare_partition_function(metric = metric, n_items = n_items)$cardinalities

# test with random sampler
n_particles <- 2
mcmc_kernel_app <- 5
num_new_obs <- 10
timesteps <- n_users / num_new_obs
timesteps2 <- dim(test_dataset)[3]
sample_dataset <- example_dataset

# run smc new user with uniform
set.seed(994)
smc_test_new_user_unif <- smc_mallows_new_users(
  rankings = sample_dataset,
  type = "partial",
  metric = metric,
  leap_size = leap_size,
  n_particles = n_particles,
  timesteps = timesteps,
  mcmc_kernel_app = mcmc_kernel_app,
  num_new_obs = num_new_obs,
  alpha_prop_sd = 0.5,
  lambda = 0.1,
  alpha_max = 20,
  aug_method = "random"
)

# run smc updated rankings with alpha unknown
smc_test_partial_unif1 <- smc_mallows_new_item_rank(
  alpha = 2,
  rankings = test_dataset, metric = metric, leap_size = leap_size,
  n_particles = n_particles,
  mcmc_kernel_app = mcmc_kernel_app, aug_method = "random",
  rho_samples_init = smc_test_new_user_unif$rho_samples[, , timesteps + 1],
  aug_rankings_init = smc_test_new_user_unif$augmented_rankings,
  alpha_fixed = TRUE
)
test_that("Updated item rank output is OK", {
  expect_is(smc_test_partial_unif1, "SMCMallows")
  expect_length(smc_test_partial_unif1, 3)
  expect_equal(dim(smc_test_partial_unif1$rho_samples), c(n_particles, n_items, 6))
  expect_length(smc_test_partial_unif1$ESS, timesteps2)
  expect_equal(dim(smc_test_partial_unif1$augmented_rankings), c(n_users, n_items, n_particles))
})

# run smc updated rankings with alpha unknown
smc_test_partial_unif2 <- smc_mallows_new_item_rank(
  rankings = test_dataset, metric = metric, leap_size = leap_size,
  n_particles = n_particles,
  mcmc_kernel_app = mcmc_kernel_app, alpha_prop_sd = 0.5,
  lambda = 0.1, alpha_max = 20, aug_method = "random",
  alpha_samples_init = smc_test_new_user_unif$alpha_samples[, timesteps + 1],
  rho_samples_init = smc_test_new_user_unif$rho_samples[, , timesteps + 1],
  aug_rankings_init = smc_test_new_user_unif$augmented_rankings
)
test_that("Updated item rank output (alpha variable) is OK", {
  expect_is(smc_test_partial_unif2, "SMCMallows")
  expect_length(smc_test_partial_unif2, 4)
  expect_equal(dim(smc_test_partial_unif2$rho_samples), c(n_particles, n_items, 6))
  expect_length(smc_test_partial_unif2$ESS, timesteps2)
  expect_equal(dim(smc_test_partial_unif2$augmented_rankings), c(n_users, n_items, n_particles))
  expect_equal(dim(smc_test_partial_unif2$alpha_samples), c(n_particles, 6))
})

# test with pseudolikelihood

smc_test_new_user_pseudo <- smc_mallows_new_users(
  rankings = example_dataset, metric = metric,
  leap_size = leap_size,
  n_particles = n_particles, timesteps = timesteps,
  mcmc_kernel_app = mcmc_kernel_app, num_new_obs = num_new_obs,
  alpha_prop_sd = 0.5, lambda = 0.1,
  alpha_max = 20, type = "partial", aug_method = "pseudolikelihood"
)

smc_test_partial_pseudo1 <- smc_mallows_new_item_rank(
  alpha = 2,
  rankings = test_dataset, metric = metric, leap_size = leap_size,
  n_particles = n_particles,
  mcmc_kernel_app = mcmc_kernel_app, aug_method = "pseudolikelihood",
  rho_samples_init = smc_test_new_user_pseudo$rho_samples[, , timesteps + 1],
  aug_rankings_init = smc_test_new_user_pseudo$augmented_rankings,
  alpha_fixed = TRUE
)
test_that("Updated item rank output is OK", {
  expect_is(smc_test_partial_pseudo1, "SMCMallows")
  expect_length(smc_test_partial_pseudo1, 3)
  expect_equal(dim(smc_test_partial_pseudo1$rho_samples), c(n_particles, n_items, 6))
  expect_length(smc_test_partial_pseudo1$ESS, timesteps2)
  expect_equal(dim(smc_test_partial_pseudo1$augmented_rankings), c(n_users, n_items, n_particles))
})

smc_test_partial_pseudo2 <- smc_mallows_new_item_rank(
  rankings = test_dataset, metric = metric, leap_size = leap_size,
  n_particles = n_particles,
  mcmc_kernel_app = mcmc_kernel_app, alpha_prop_sd = 0.5,
  lambda = 0.1, alpha_max = 20, aug_method = "pseudolikelihood",
  alpha_samples_init = smc_test_new_user_unif$alpha_samples[, timesteps + 1],
  rho_samples_init = smc_test_new_user_unif$rho_samples[, , timesteps + 1],
  aug_rankings_init = smc_test_new_user_unif$augmented_rankings
)
test_that("Updated item rank output (variable alpha) is OK", {
  expect_is(smc_test_partial_pseudo2, "SMCMallows")
  expect_length(smc_test_partial_pseudo2, 4)
  expect_equal(dim(smc_test_partial_pseudo2$rho_samples), c(n_particles, n_items, 6))
  expect_length(smc_test_partial_pseudo2$ESS, timesteps2)
  expect_equal(dim(smc_test_partial_pseudo2$augmented_rankings), c(n_users, n_items, n_particles))
  expect_equal(dim(smc_test_partial_pseudo2$alpha_samples), c(n_particles, 6))
})

# check metric and aug_method error
test_that("metric and aug_method must match", {
  expect_error(
    smc_mallows_new_item_rank(
      alpha = 2,
      rankings = test_dataset, metric = "cayley", leap_size = leap_size,
      n_particles = n_particles,
      mcmc_kernel_app = mcmc_kernel_app, aug_method = "pseudolikelihood",
      alpha_fixed = TRUE
    ),
    "Pseudolikelihood only supports footrule and spearman metrics"
  )
})
