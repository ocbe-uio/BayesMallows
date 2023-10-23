# Assume our first dataset contains 1000 preferences, observed in 10 batches
data_batch1 <- sushi_rankings[1:1000, ]
num_new_obs <- 100

# Calculating rho and alpha samples
samples <- smc_mallows_new_users(
  rankings = data,
  n_particles = 1000,
  timesteps = 10,
  mcmc_kernel_app = 10,
  num_new_obs = num_new_obs,
  verbose = TRUE
)

# Studying the structure of the output
plot(samples)


