# Reproduce Section 6.3.2 of Anja's PhD thesis
library(BayesMallows)
dat <- sushi_rankings[1:100, ]

# Initialize
dat_init <- dat
dat_init[dat_init > 5] <- NA

mod_init <- smc_mallows_new_users(
  rankings = dat_init,
  type = "partial",
  n_particles = 1000,
  timesteps = 10,
  mcmc_kernel_app = 10,
  num_new_obs = 10,
  aug_method = "pseudolikelihood"
)

# Add particles


rankings <- array(dim = c(dim(dat), 40))

assessor_group <- rep(1:10, each = 10)
inds <- expand.grid(group = 1:10, k = 6:9)

for(i in seq_len(nrow(inds))) {
  k <- inds[i, "k"]
  group <- inds[i, "group"]
  rankings[,, i] <- dat
  rankings[assessor_group <= i,, i][rankings[assessor_group <= i,, i] < k] <- NA
}

mod_item <- smc_mallows_new_item_rank(
  rankings = rankings[,,1:20],
  n_particles = 1000,
  mcmc_kernel_app = 10,
  aug_rankings_init = mod_init$augmented_rankings,
  rho_samples_init = mod_init$rho_samples[,, 11],
  alpha_samples_init = mod_init$alpha_samples[, 11],
  aug_method = "pseudolikelihood",
  verbose = TRUE
)

apply(mod_item$alpha_samples, 2, mean)
