# Assume our first dataset contains 100 preferences, observed in 10 batches
data_batch1 <- sushi_rankings[1:100, ]

# Calculating rho and alpha samples
mod1 <- smc_mallows_new_users(
  rankings = data_batch1,
  n_particles = 1000,
  timesteps = 1,
  mcmc_kernel_app = 10,
  num_new_obs = 100
)

# We can plot the posterior of alpha
plot(mod1)

# Next, assume a new batch arrives, with 100 new preferences
data_batch2 <- sushi_rankings[101:200, ]

# We update the model in a single timestep. By default, the settings from the
# call to smc_mallows_new_users() is retained.
mod2 <- smc_mallows_update(
  model = mod1,
  rankings = data_batch2
)

plot(mod2)
