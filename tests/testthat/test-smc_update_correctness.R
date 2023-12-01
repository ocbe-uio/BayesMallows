test_that("update_mallows is correct for new rankings", {
  triple_potato <- rbind(potato_visual, potato_visual, potato_visual)
  rownames(triple_potato) <- seq_len(nrow(triple_potato))
  user_ids <- rownames(triple_potato)

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(triple_potato),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  mod_init <- compute_mallows(
    data = setup_rank_data(triple_potato[1:4, , drop = FALSE]),
    compute_options = set_compute_options(nmc = 20000, burnin = 1000)
  )

  mod_smc <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = triple_potato[5:20, ]),
    smc_options = set_smc_options(n_particles = 10000, mcmc_steps = 15)
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(rankings = triple_potato[21:36, ])
  )

  expect_equal(mean(mod_smc_next$alpha$value),
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    tolerance = 0.01
  )

  expect_equal(sd(mod_smc_next$alpha$value),
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    tolerance = 0.1
  )

  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  expect_lte(
    compute_rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    1
  )

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(sushi_rankings),
    compute_options = set_compute_options(nmc = 2000, burnin = 200)
  )

  mod_init <- compute_mallows(
    data = setup_rank_data(sushi_rankings[1:100, ]),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  mod_smc <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = sushi_rankings[101:2000, ]),
    smc_options = set_smc_options(n_particles = 500, mcmc_steps = 5)
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(sushi_rankings[2001:5000, ])
  )

  # Posterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_smc_next$alpha$value),
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
    tolerance = 0.02
  )

  expect_equal(sd(mod_smc_next$alpha$value),
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
    tolerance = 0.05
  )

  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  expect_lte(
    compute_rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    2
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

  for (aug in c("uniform", "pseudo")) {
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

    expect_equal(
      mean(mod_smc_next$alpha$value),
      mean(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 1000]),
      tolerance = .02
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
      ), 3
    )
  }
})

test_that("update_mallows is correct for updated partial rankings", {
  set.seed(123)
  user_ids <- rownames(potato_visual)
  dat0 <- potato_visual
  dat0[] <- ifelse(runif(length(dat0)) > .5, NA_real_, dat0)

  mod0 <- compute_mallows(
    data = setup_rank_data(rankings = dat0),
    compute_options = set_compute_options(nmc = 10000, burnin = 5000)
  )

  dat1 <- potato_visual
  dat1 <- ifelse(is.na(dat0) & runif(length(dat1)) > .5, NA_real_, dat1)

  mod1 <- update_mallows(
    model = mod0,
    new_data = setup_rank_data(rankings = dat1, user_ids = user_ids),
    smc_options = set_smc_options(n_particles = 10000, mcmc_steps = 10)
  )

  mod_bmm1 <- compute_mallows(
    data = setup_rank_data(rankings = dat1),
    compute_options = set_compute_options(nmc = 50000)
  )

  expect_equal(
    mean(mod1$alpha$value),
    mean(mod_bmm1$alpha$value[mod_bmm1$alpha$iteration > 5000]),
    tolerance = .1
  )

  dat2 <- potato_visual
  dat2 <- ifelse(is.na(dat1) & runif(length(dat2)) > .5, NA_real_, dat2)

  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(rankings = dat2, user_ids = user_ids)
  )

  mod_bmm <- compute_mallows(
    data = setup_rank_data(rankings = dat2),
    compute_options = set_compute_options(nmc = 50000)
  )

  expect_equal(
    mean(mod2$alpha$value),
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 5000]),
    tolerance = .1
  )

  expect_equal(
    sd(mod2$alpha$value),
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 5000]),
    tolerance = .1
  )
})

test_that("update_mallows is correct for new top-k rankings", {
  skip_on_cran()

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
    smc_options = set_smc_options(n_particles = 10000, mcmc_steps = 10)
  )

  mod_bmm1 <- compute_mallows(
    data = setup_rank_data(rankings = dat1),
    compute_options = set_compute_options(nmc = 50000)
  )

  expect_equal(
    mean(mod1$alpha$value),
    mean(mod_bmm1$alpha$value[mod_bmm1$alpha$iteration > 5000]),
    tolerance = .1
  )

  expect_equal(
    sd(mod1$alpha$value),
    sd(mod_bmm1$alpha$value[mod_bmm1$alpha$iteration > 5000]),
    tolerance = .1
  )

  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(rankings = dat2, user_ids = user_ids)
  )

  mod_bmm2 <- compute_mallows(
    data = setup_rank_data(rankings = dat2),
    compute_options = set_compute_options(nmc = 50000)
  )

  expect_equal(
    mean(mod2$alpha$value),
    mean(mod_bmm2$alpha$value[mod_bmm2$alpha$iteration > 5000]),
    tolerance = .1
  )

  expect_equal(
    sd(mod2$alpha$value),
    sd(mod_bmm2$alpha$value[mod_bmm2$alpha$iteration > 5000]),
    tolerance = .1
  )
})
