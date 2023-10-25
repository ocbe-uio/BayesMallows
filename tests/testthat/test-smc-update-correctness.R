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

test_that("smc_mallows_new_users is correct for new partial rankings", {
  skip_extended()
  set.seed(123)

  rankings <- matrix(rep(c(
    1, 2, 3,
    1, 3, 2,
    1, 2, 3,
    1, 2, 3,
    2, 1, 3), times = 10), ncol = 3, byrow = TRUE)

  rankings[sample(seq_along(rankings), 10)] <- NA

  bmm_mod <- compute_mallows(rankings = rankings)
  bmm_mod$burnin <- 100

  smc_onego <- smc_mallows_new_users(
    rankings = rankings,
    type = "partial",
    n_particles = 1000,
    timesteps = 10,
    mcmc_kernel_app = 10,
    num_new_obs = 5,
    verbose = TRUE
  )

  inds <- rep(1:10, each = 5)
  smc_init <- smc_mallows_new_users(
    rankings = rankings[inds == 1, ],
    type = "partial",
    n_particles = 1000,
    timesteps = 1,
    mcmc_kernel_app = 10,
    num_new_obs = 5,
    verbose = TRUE
  )

  smc_update <- smc_init
  for(i in 2:10) {
    smc_update <- smc_mallows_update(
      model = smc_update, rankings = rankings[inds == i, ],
      verbose = TRUE
    )
  }

  expect_equal(mean(smc_update$alpha_samples[, 2]), 2.51, tolerance = .01)
  expect_equal(mean(smc_onego$alpha_samples[, 11]), 2.50, tolerance = .01)
  expect_equal(mean(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 100]), 2.48,
               tolerance = .01)

  expect_equal(sd(smc_update$alpha_samples[, 2]), .34, tolerance = .01)
  expect_equal(sd(smc_onego$alpha_samples[, 11]), .35, tolerance = .01)
  expect_equal(sd(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 100]), 0.34,
               tolerance = .01)

  expect_equal(consensus_unwrapper(smc_update), 1:3)
  expect_equal(consensus_unwrapper(smc_onego), 1:3)
  expect_equal(consensus_unwrapper(bmm_mod), 1:3)

})

