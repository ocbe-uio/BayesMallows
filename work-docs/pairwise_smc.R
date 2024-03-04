library(BayesMallows)
library(patchwork)
set.seed(1)

dat <- subset(beach_preferences, assessor == 1)

mod <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 3000)
)

for(i in 2:60) {
  print(i)
  dat <- subset(beach_preferences, assessor == i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(preferences = dat, timepoint = i),
    smc_options = set_smc_options(latent_sampling_lag = 0, mcmc_steps = 20)
  )
}

mod_bmm <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 10000, burnin = 5000)
)

mod_smc <- update_mallows(
  model = mod_bmm,
  new_data = setup_rank_data(preferences = beach_preferences),
  smc_options = set_smc_options(n_particles = 20, mcmc_steps = 3000)
)

plot(mod_smc)

plot(mod) + plot(mod_bmm)

plot(mod, parameter = "rho", items = 1:15) +
  plot(mod_bmm, parameter = "rho", items = 1:15)

