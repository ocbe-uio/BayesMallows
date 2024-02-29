devtools::load_all()

dat0 <- subset(beach_preferences, assessor %in% 1:3)

mod0 <- compute_mallows(
  data = setup_rank_data(preferences = dat0),
  compute_options = set_compute_options(nmc = 10000, burnin = 3000)
)
assess_convergence(mod0)

dat1 <- subset(beach_preferences, assessor %in% 4:6)
mod1 <- update_mallows(
  model = mod0,
  new_data = setup_rank_data(preferences = dat1),
  smc_options = set_smc_options(n_particles = 2000, mcmc_steps = 5)
)

plot(mod0)
plot(mod1)
