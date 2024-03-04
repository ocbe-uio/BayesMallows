library(BayesMallows)
library(patchwork)
set.seed(1)

dat <- subset(beach_preferences, assessor < 60)

mod_init <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

assess_convergence(mod_init)
dat_60 <- subset(beach_preferences, assessor == 60)

mod <- update_mallows(
  model = mod_init,
  new_data = setup_rank_data(
    preferences = dat_60[1, , drop = FALSE], n_items = 15),
  smc_options = set_smc_options(n_particles = 4, mcmc_steps = 2)
)

plot(mod)
