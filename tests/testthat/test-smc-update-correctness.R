
test_that("skip_extended works", {
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "true")
  expect_condition(skip_extended(), NA, class = "skip")
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "")
  expect_condition(skip_extended(), class = "skip")
})


test_that("smc_mallows_update is correct", {
  skip_extended()
  set.seed(123)
  mod_ref <- compute_mallows(potato_visual, nmc = 20000)
  mod_ref$burnin <- 1000

  mod_smc <- smc_mallows_new_users(
    potato_visual[1, , drop = FALSE],
    timesteps = 1,
    num_new_obs = 1,
    n_particles = 5000,
    mcmc_kernel_app = 10
  )

  for(i in seq(from = 2, to = nrow(potato_visual))) {
    mod_smc <- smc_mallows_update(mod_smc, potato_visual[i, , drop = FALSE])
  }

  c1 <- compute_consensus(mod_ref)
  c2 <- compute_consensus(mod_smc)

  # Is there any disagreement between the methods about the ranking of the items?
  expect_equal(
    grep("[0-9]+", c1$item),
    grep("[0-9]+", c2$item)
  )

})
