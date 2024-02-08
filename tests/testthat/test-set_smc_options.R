test_that("set_smc_options works", {
  expect_error(
    set_smc_options(n_particles = 3.4),
    "n_particles must be a positive integer")
  expect_error(
    set_smc_options(n_particles = -1),
    "n_particles must be a positive integer")
  expect_error(
    set_smc_options(mcmc_steps = 3.4),
    "mcmc_steps must be a positive integer")
  expect_error(
    set_smc_options(mcmc_steps = -1),
    "mcmc_steps must be a positive integer")
  expect_error(
    set_smc_options(latent_sampling_lag = -1),
    "latent_sampling_lag must be a positive integer"
  )

  val <- set_smc_options(n_particles = 35, mcmc_steps = 2)
  expect_equal(val$n_particles, 35)
  expect_equal(val$mcmc_steps, 2)
  expect_equal(val$latent_sampling_lag, NA_integer_)

  val <- set_smc_options(latent_sampling_lag = 3)
  expect_equal(val$n_particles, 1000)
  expect_equal(val$mcmc_steps, 5)
  expect_equal(val$latent_sampling_lag, 3)
})
