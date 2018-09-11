# Under construction
# Number of items
n_items <- 15
# Set the consensus ranking
rho0 <- seq(from = 1, to = n_items, by = 1)
# Set the scale
alpha0 <- 5
# Number of samples we want
n_samples <- 2
# Burnin
burnin <- 1000
# Thinning
thinning <- 100

sample_mallows(rho0, alpha0, n_samples, burnin, thinning,
               diagnostic = TRUE)
