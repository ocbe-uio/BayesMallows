# Example dataset potato_partial has 11 timepoints with incrementally updated
# rankings
updated_partial_mod <- smc_mallows_new_item_rank(
  rankings = potato_partial,
  n_particles = 100,
  mcmc_kernel_app = 10,
  aug_method = "pseudolikelihood"
)

# Plot posterior of alpha
plot(updated_partial_mod)
