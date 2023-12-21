library(BayesMallows)

set.seed(123)
dat0 <- t(apply(potato_visual, 1, function(x) {
  inds <- sample(length(x), 2)
  x[inds] <- NA
  x
}))
dat <- rbind(dat0, dat0, dat0)
rownames(dat) <- seq_len(nrow(dat))
user_ids <- rownames(dat)

bmm_mod <- compute_mallows(
  data = setup_rank_data(dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

mod_init <- compute_mallows(
  data = setup_rank_data(dat[1:4, ]),
  compute_options = set_compute_options(
    nmc = 10000, burnin = 1000,
    save_aug = TRUE
  )
)

aug <- "pseudo"
system.time({
  mod_smc <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = dat[5:20, ]),
    smc_options = set_smc_options(
      n_particles = 10000, mcmc_steps = 10, aug_method = aug
    )
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(rankings = dat[21:36, ])
  )
})
# 4.7 seconds


set.seed(1)
beach_data <- setup_rank_data(
  preferences = beach_preferences
)
system.time({
  model_fit <- compute_mallows(
    data = beach_data,
    compute_options = set_compute_options(save_aug = TRUE))
})
# 0.447 seconds
# 0.33 seconds
