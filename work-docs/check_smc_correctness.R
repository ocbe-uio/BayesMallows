library(BayesMallows)

mod_ref <- compute_mallows(potato_weighing)

mod_smc <- smc_mallows_new_users(
  potato_weighing[1, , drop = FALSE],
  timesteps = 1,
  num_new_obs = 1,
  n_particles = 100,
  mcmc_kernel_app = 5
  )

for(i in seq(from = 2, to = nrow(potato_weighing))) {
  mod_smc <- smc_mallows_update(mod_smc, potato_weighing[i, , drop = FALSE])
}

plot(mod_smc)
plot(mod_ref, burnin = 500)

