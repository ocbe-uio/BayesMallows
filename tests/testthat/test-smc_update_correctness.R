test_that("update_mallows is correct for new rankings", {

  triple_potato <- rbind(potato_visual, potato_visual, potato_visual)

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(triple_potato),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000))

  mod_init <- compute_mallows(
    data = setup_rank_data(triple_potato[1:4, , drop = FALSE]),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000))

  mod_smc <- update_mallows(
    model = mod_init,
    new_rankings = triple_potato[5:20, ],
    n_particles = 5000,
    mcmc_steps = 3
    )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_rankings = triple_potato[21:36, ]
  )


  # Poterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_smc_next$alpha$value),
               mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
               tolerance = 0.01)

  expect_equal(sd(mod_smc_next$alpha$value),
               sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
               tolerance = 0.1)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  # No items should be in disagreement
  expect_equal(
    compute_rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    0
  )

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(sushi_rankings),
    compute_options = set_compute_options(nmc = 1000, burnin = 200))

  mod_init <- compute_mallows(
    data = setup_rank_data(sushi_rankings[1:100, ]),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000))

  mod_smc <- update_mallows(
    model = mod_init,
    new_rankings = sushi_rankings[101:2000, ],
    n_particles = 500,
    mcmc_steps = 3
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_rankings = sushi_rankings[2001:5000, ]
  )

  # Poterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_smc_next$alpha$value),
               mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
               tolerance = 0.01)

  expect_equal(sd(mod_smc_next$alpha$value),
               sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
               tolerance = 0.1)

  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  # How many items are in disagreement
  expect_equal(
    compute_rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    0
  )

})

test_that("update_mallows is correct for new partial rankings", {
  skip_on_cran()

  set.seed(123)

  dat0 <- t(apply(potato_visual, 1, function(x) {
    inds <- sample(length(x), 2)
    x[inds] <- NA
    x
  }))
  dat <- rbind(dat0, dat0, dat0)

  bmm_mod <- compute_mallows(
    data = setup_rank_data(dat),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  mod_init <- compute_mallows(
    data = setup_rank_data(dat[1:4, ]),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000,
                                          save_aug = TRUE))

  for(aug in c("uniform", "pseudo")) {
    mod_smc <- update_mallows(
      model = mod_init,
      new_rankings = dat[5:20, ],
      n_particles = 10000,
      mcmc_steps = 10,
      augmentation = aug
    )

    mod_smc_next <- update_mallows(
      model = mod_smc,
      new_rankings = dat[21:36, ]
    )

    expect_equal(
      mean(mod_smc_next$alpha$value),
      mean(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 1000]),
      tolerance = ifelse(aug == "uniform", .2, .05)
    )

    expect_equal(
      sd(mod_smc_next$alpha$value),
      sd(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 1000]),
      tolerance = .2
    )

    bmm_consensus <- compute_consensus(bmm_mod)
    smc_consensus <- compute_consensus(mod_smc_next)

    expect_lte(
      compute_rank_distance(
        matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
        as.numeric(as.factor(smc_consensus$item)),
        metric = "ulam"
      ),
      ifelse(aug == "uniform", 3, 1)
    )
  }

})
