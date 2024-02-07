
library(patchwork)
devtools::load_all()
dat <- potato_visual
#dat[dat > 15] <- NA

mod_bmm <- compute_mallows(
  data = setup_rank_data(dat),
  compute_options = set_compute_options(nmc = 10000))
assess_convergence(mod_bmm)

mod_t <- sample_prior(10000, 20, set_priors(lambda = .1))

for(i in 1:12) {
  print(i)
  mod_t <- update_mallows(
    model = mod_t,
    new_data = setup_rank_data(rankings = dat[i, ]),
    smc_options = set_smc_options(n_particles = 10000, mcmc_steps = 10, latent_sampling_lag = 1)
  )
}


plot(mod_t) + plot(mod_bmm, burnin = 5000)
