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
    preferences = dat_60[1, , drop = FALSE], n_items = 15, user_ids = 60)
)


for(i in seq(from = 2, to = nrow(dat_60))) {
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(
      preferences = dat_60[1:i, , drop = FALSE], n_items = 15, user_ids = 60
    )
  )
}

plot(mod)

mod <- update_mallows(
  model = mod_init,
  new_data = setup_rank_data(preferences = beach_preferences),
  smc_options = set_smc_options(n_particles = 20, mcmc_steps = 10000)
)
plot(mod)
