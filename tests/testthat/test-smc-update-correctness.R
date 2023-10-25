consensus_unwrapper <- function(x) {
  x <- compute_consensus(x)
  as.numeric(unlist(regmatches(x$item, gregexpr("[0-9]+", x$item))))
}


test_that("skip_extended works", {
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "true")
  expect_condition(skip_extended(), NA, class = "skip")
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "")
  expect_condition(skip_extended(), class = "skip")
})


test_that("smc_mallows_update is correct for new rankings", {
  skip_extended()
  set.seed(123)
  # Metropolis-Hastings
  mod_bmm <- compute_mallows(rankings = potato_visual, nmc = 10000)
  mod_bmm$burnin <- 1000

  # Single function call, taking one observation at a time
  mod_ref <- smc_mallows_new_users(
    rankings = potato_visual,
    n_particles = 1000,
    timesteps = nrow(potato_visual),
    num_new_obs = 1,
    mcmc_kernel_app = 50
  )

  # Sequentially, using update function
  mod_smc <- smc_mallows_new_users(
    potato_visual[1, , drop = FALSE],
    timesteps = 1,
    num_new_obs = 1,
    n_particles = 1000,
    mcmc_kernel_app = 50
  )

  for(i in seq(from = 2, to = nrow(potato_visual))) {
    mod_smc <- smc_mallows_update(model = mod_smc,
                                  rankings = potato_visual[i, , drop = FALSE])
  }

  # Poterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_smc$alpha_samples[, 2]), 10.8, tolerance = 1e-2)
  expect_equal(mean(mod_ref$alpha_samples[, 13]), 10.8, tolerance = 1e-2)
  expect_equal(mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]), 10.86, tolerance = 1e-2)

  # Same test for posterior standard deviation
  expect_equal(sd(mod_smc$alpha_samples[, 2]), 0.73, tolerance = 1e-2)
  expect_equal(sd(mod_ref$alpha_samples[, 13]), 0.73, tolerance = 1e-2)
  expect_equal(sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]), 0.77, tolerance = 1e-2)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- consensus_unwrapper(mod_bmm)
  ref_consensus <- consensus_unwrapper(mod_ref)
  smc_consensus <- consensus_unwrapper(mod_smc)

  # How many items are in disagreement
  expect_equal(
    rank_distance(matrix(ref_consensus, nrow = 1),
                  bmm_consensus, metric = "ulam"),
    0)
  expect_equal(
    rank_distance(matrix(ref_consensus, nrow = 1),
                  smc_consensus, metric = "ulam"),
    1)



})

test_that("smc_mallows_update is correct for new partial rankings", {
  # Introduce missingness in potato_visual
  set.seed(123)
  dat <- potato_visual
  dat[dat > 15] <- NA

  # Metropolis-Hastings
  mod_bmm <- compute_mallows(rankings = dat, nmc = 10000)
  mod_bmm$burnin <- 1000

  # Single function call, taking one observation at a time
  mod_ref <- smc_mallows_new_users(
    rankings = dat,
    type = "partial",
    n_particles = 2000,
    timesteps = nrow(dat),
    num_new_obs = 1,
    mcmc_kernel_app = 50,
    aug_method = "pseudolikelihood"
  )
  mod_ref_uniform <- smc_mallows_new_users(
    rankings = dat,
    type = "partial",
    n_particles = 2000,
    timesteps = nrow(dat),
    num_new_obs = 1,
    mcmc_kernel_app = 50,
    aug_method = "random"
  )

  # Poterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_ref$alpha_samples[, 13]), 10.6, tolerance = 0.01)
  expect_equal(mean(mod_ref_uniform$alpha_samples[, 13]), 10.4, tolerance = 0.01)
  expect_equal(mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]), 10.8, tolerance = 0.01)

  # Same test for posterior standard deviation
  expect_equal(sd(mod_ref$alpha_samples[, 13]), 0.73, tolerance = 0.1)
  expect_equal(sd(mod_ref_uniform$alpha_samples[, 13]), 0.74, tolerance = 0.1)
  expect_equal(sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]), 0.77, tolerance = 0.1)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- consensus_unwrapper(mod_bmm)
  ref_consensus <- consensus_unwrapper(mod_ref)
  ref_uniform_consensus <- consensus_unwrapper(mod_ref_uniform)

  # How many items are in disagreement
  expect_equal(
    rank_distance(matrix(ref_consensus, nrow = 1),
                  bmm_consensus, metric = "ulam"),
    5)
  expect_equal(
    rank_distance(matrix(ref_uniform_consensus, nrow = 1),
                  bmm_consensus, metric = "ulam"),
    3)

})

