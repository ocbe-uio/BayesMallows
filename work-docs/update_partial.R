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
  mcmc_steps = 10,
  num_new_obs = 5,
  verbose = TRUE
)

inds <- rep(1:10, each = 5)
smc_init <- smc_mallows_new_users(
  rankings = rankings[inds == 1, ],
  type = "partial",
  n_particles = 1000,
  timesteps = 1,
  mcmc_steps = 10,
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

