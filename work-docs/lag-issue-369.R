rm(list=ls())
library(patchwork)
devtools::load_all()
dat <- potato_visual
dat[dat < 15] <- NA

mod0 <- compute_mallows(
  data = setup_rank_data(dat[1:3, ]),
  compute_options = set_compute_options(nmc = 10000, burnin = 100)
  )

timepoint <- 1
mod1 <- mod0
for(i in 4:12) {
  mod1 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(rankings = dat[i, ]),
    smc_options = set_smc_options(n_particles = 2000)
  )
  timepoint <- timepoint + 1
}


timepoint <- 1
mod2 <- mod0
for(i in 4:12) {
  mod2 <- update_mallows(
    model = mod2,
    new_data = setup_rank_data(rankings = dat[i, ], timepoint = timepoint),
    smc_options = set_smc_options(n_particles = 2000, latent_sampling_lag = 1)
  )
  timepoint <- timepoint + 1
}

plot(mod1) + plot(mod2)
