library(BayesMallows)
sessionInfo()
library(ggplot2)
dat0 <- ifelse(potato_visual > 10, NA_real_, potato_visual)
dat1 <- ifelse(potato_visual > 12, NA_real_, potato_visual)
dat2 <- ifelse(potato_visual > 14, NA_real_, potato_visual)

set.seed(123)
user_ids <- rownames(potato_visual)

mod0 <- compute_mallows(
  data = setup_rank_data(rankings = dat0),
  compute_options = set_compute_options(nmc = 10000, burnin = 5000)
)

mod1 <- update_mallows(
  model = mod0,
  new_data = setup_rank_data(rankings = dat1, user_ids = user_ids),
  compute_options = set_compute_options(leap_size = 4),
  smc_options = set_smc_options(n_particles = 10000, mcmc_steps = 10)
)

plot(mod1) + ggtitle("SMC posterior with dat1")

mod_bmm1 <- compute_mallows(
  data = setup_rank_data(rankings = dat1),
  compute_options = set_compute_options(nmc = 50000)
)

plot(mod_bmm1, burnin = 5000) + ggtitle("MCMC posterior with dat1")

mod2 <- update_mallows(
  model = mod1,
  new_data = setup_rank_data(rankings = dat2, user_ids = user_ids)
)

plot(mod2) + ggtitle("SMC posterior with dat2")

mod_bmm <- compute_mallows(
  data = setup_rank_data(rankings = dat2),
  compute_options = set_compute_options(nmc = 50000)
)

plot(mod_bmm, burnin = 5000) + ggtitle("MCMC posterior with dat2")
