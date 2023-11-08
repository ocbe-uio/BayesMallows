set.seed(123)
mod_ref <- compute_mallows(potato_visual, nmc = 20000)
mod_ref$burnin <- 1000

mod_smc <- smc_mallows_new_users(
  potato_visual[1, , drop = FALSE],
  timesteps = 1,
  num_new_obs = 1,
  n_particles = 5000,
  mcmc_steps = 10
  )

for(i in seq(from = 2, to = nrow(potato_visual))) {
  mod_smc <- smc_mallows_update(mod_smc, potato_visual[i, , drop = FALSE])
}

plot(mod_smc)
plot(mod_ref, burnin = 1000)

mean(mod_smc$alpha_samples[, 2])
mean(mod_ref$alpha$value[mod_ref$alpha$iteration > 1000])

c1 <- compute_consensus(mod_ref)
c2 <- compute_consensus(mod_smc)

# Is there any disagreement between the methods about the ranking of the items?
mean(grep("[0-9]+", c1$item) != grep("[0-9]+", c2$item))
