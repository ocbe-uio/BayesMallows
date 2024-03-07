library(BayesMallows)
library(patchwork)
set.seed(1)

dat <- subset(beach_preferences, assessor < 5)

mod_init <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)
# assess_convergence(mod_init)
alpha <- numeric()

mod <- mod_init
for(i in 6:60) {
  print(i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(preferences = subset(beach_preferences, assessor == i),
                               timepoint = i),
    smc_options = set_smc_options(
      n_particles = 10000, mcmc_steps = 100)
  )
  alpha <- c(alpha, mean(mod$alpha_samples))
}

#plot(alpha)

mod_bmm <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

plot(mod_bmm) + plot(mod) + plot_layout(ncol = 1)
