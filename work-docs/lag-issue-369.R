rm(list=ls())
devtools::load_all()
dat <- potato_visual
dat[dat < 15] <- NA

mod0 <- compute_mallows(
  data = setup_rank_data(dat[1:2, ]),
  compute_options = set_compute_options(burnin = 500)
  )

mod1 <- update_mallows(
  model = mod0,
  new_data = setup_rank_data(rankings = dat[3:4, ], timepoint = c(1, 1)),
  smc_options = set_smc_options(latent_sampling_lag = 20)
  )

mod2 <- update_mallows(
  model = mod1,
  new_data = setup_rank_data(rankings = dat[5:6, ], timepoint = c(2, 2))
)
