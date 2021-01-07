# Create example data. We set the burn-in and thinning very low
# for the sampling to go fast
data0 <- sample_mallows(rho0 = 1:5, alpha=10, n_samples = 1000,
                        burnin = 10, thinning = 1)
# Find the frequency distribution
rank_freq_distr(rankings=data0)
