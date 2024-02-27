rm(list=ls())
devtools::load_all()
ids <- rownames(potato_visual)

mod0 <- sample_prior(100, 20)

dat1 <- potato_visual
dat1[dat1 > 18] <- NA
part1 <- 1:2
mod1 <- update_mallows(
  model = mod0,
  new_data = setup_rank_data(rankings = dat1[ids[part1], ], user_ids = ids[part1]),
  smc_options = set_smc_options(n_particles = 2)
  )

dat2 <- potato_visual
dat2[dat2 > 19] <- NA
part2 <- 2:3

mod2 <- update_mallows(
  model = mod1,
  new_data = setup_rank_data(rankings = dat2[ids[part2], ], user_ids = ids[part2])
)

