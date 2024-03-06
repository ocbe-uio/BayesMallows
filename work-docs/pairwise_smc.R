library(BayesMallows)
library(patchwork)
set.seed(1)

dat <- subset(beach_preferences, assessor < 5)

mod_init <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)
assess_convergence(mod_init)
alpha <- numeric()

nd <- setup_rank_data(
  preferences = beach_preferences[1:3, , drop = FALSE], n_items = 15
)

mod <- mod_init
for(i in 6:20) {
  print(i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(preferences = subset(beach_preferences, assessor == i),
                               timepoint = i),
    smc_options = set_smc_options(latent_sampling_lag = 0)
  )
  alpha <- c(alpha, mean(mod$alpha_samples))
}

plot(alpha)

assess_convergence(mod_init)
dat_60 <- subset(beach_preferences, assessor == 60)

nd <- setup_rank_data(
  preferences = dat_60[1:5, , drop = FALSE],
  n_items = 15,
  user_ids = 60,
  shuffle_unranked = TRUE
)
nd$rankings

mod <- update_mallows(
  model = mod_init,
  new_data = setup_rank_data(
    preferences = dat_60[, , drop = FALSE], n_items = 15, user_ids = 60, shuffle_unranked = TRUE),
  smc_options = set_smc_options(mcmc_steps = 15)
)

plot(mod) + plot(mod_init)

for(i in seq(from = 2, to = 10, by = 1)) {
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
