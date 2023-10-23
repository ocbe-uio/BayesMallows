
test_that("skip_extended works", {
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "true")
  expect_condition(skip_extended(), NA, class = "skip")
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "")
  expect_condition(skip_extended(), class = "skip")
})


test_that("smc_mallows_update is correct", {
  skip_extended()
  set.seed(123)
  mod_ref <- compute_mallows(potato_visual, nmc = 5000)
  mod_ref$burnin <- 1000

  mod_smc <- smc_mallows_new_users(
    potato_visual[1, , drop = FALSE],
    timesteps = 1,
    num_new_obs = 1,
    n_particles = 10000,
    mcmc_kernel_app = 10
  )

  for(i in seq(from = 2, to = nrow(potato_visual))) {
    mod_smc <- smc_mallows_update(mod_smc, potato_visual[i, , drop = FALSE])
  }

  c1 <- compute_consensus(mod_ref)
  c2 <- compute_consensus(mod_smc)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- unlist(regmatches(c1$item, gregexpr("[0-9]+", c1$item)))
  smc_consensus <- unlist(regmatches(c2$item, gregexpr("[0-9]+", c2$item)))
  expect_equal(
    bmm_consensus[1:3], smc_consensus[1:3]
  )

  # How many ranks are in disagreement
  # The results, 6, is much less than expected by chance
  expect_equal(
    rank_distance(matrix(as.numeric(bmm_consensus), nrow = 1),
                  as.numeric(smc_consensus),
                  metric = "ulam"), 6)



})
