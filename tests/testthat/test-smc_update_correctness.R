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
    new_data = setup_rank_data(rankings = triple_potato[5:20, ])
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(rankings = triple_potato[21:36, ])
  )

  expect_equal(
    mean(mod_smc_next$alpha$value),
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    tolerance = 0.01
  )

  expect_equal(
    sd(mod_smc_next$alpha$value),
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
    2
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
    tolerance = 0.01
  )

  expect_equal(sd(mod_smc_next$alpha$value),
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
    tolerance = 0.04
  )

  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  expect_lte(
    compute_rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    0
  )
})

test_that("update_mallows is correct for new partial rankings", {
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
      smc_options = set_smc_options(n_particles = 2000),
      compute_options = set_compute_options(aug_method = aug)
    )

    mod_smc_next <- update_mallows(
      model = mod_smc,
      new_data = setup_rank_data(rankings = dat[21:36, ])
    )

    expect_equal(
      mean(mod_smc_next$alpha$value),
      mean(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 1000]),
      tolerance = .01
    )

    expect_equal(
      sd(mod_smc_next$alpha$value),
      sd(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 1000]),
      tolerance = ifelse(aug == "uniform", .2, .1)
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
  set.seed(1)
  user_ids <- rownames(potato_visual)
  dat0 <- potato_visual
  dat0[] <- ifelse(runif(length(dat0)) > .5, NA_real_, dat0)

  mod0 <- compute_mallows(
    data = setup_rank_data(rankings = dat0),
    compute_options = set_compute_options(nmc = 20000, burnin = 2500)
  )

  dat1 <- potato_visual
  dat1 <- ifelse(is.na(dat0) & runif(length(dat1)) > .5, NA_real_, dat1)

  mod1 <- update_mallows(
    model = mod0,
    new_data = setup_rank_data(rankings = dat1, user_ids = user_ids),
    compute_options = set_compute_options(aug_method = "pseudo")
  )

  mod_bmm1 <- compute_mallows(
    data = setup_rank_data(rankings = dat1),
    compute_options = set_compute_options(
      nmc = 5000, burnin = 2000,
      aug_method = "pseudo"
    )
  )

  expect_equal(
    mean(mod1$alpha$value),
    mean(mod_bmm1$alpha$value[mod_bmm1$alpha$iteration > 4000]),
    tolerance = .05
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
})

test_that("update_mallows does not suffer from numerical overflow", {
  data_batch1 <- sushi_rankings[1:100, ]
  mod1 <- compute_mallows(
    data = setup_rank_data(data_batch1),
    compute_options = set_compute_options(burnin = 200)
  )

  data_batch2 <- sushi_rankings[101:200, ]
  data_batch2[data_batch2 > 3] <- NA
  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(data_batch2, user_ids = 1:100),
    smc_options = set_smc_options(n_particles = 50),
    compute_options = set_compute_options(aug_method = "pseudo")
  )

  expect_s3_class(mod2, "SMCMallows")

  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(data_batch2, user_ids = 1:100),
    smc_options = set_smc_options(n_particles = 50),
    compute_options = set_compute_options(aug_method = "uniform")
  )

  expect_s3_class(mod2, "SMCMallows")
})

test_that("update_mallows works for data one at a time", {
  n <- 50
  set.seed(2)
  mod_bmm <- compute_mallows(data = setup_rank_data(sushi_rankings[1:n, ]))

  mod <- sample_prior(1000, ncol(sushi_rankings))
  for (i in seq_len(n)) {
    mod <- update_mallows(
      model = mod,
      new_data = setup_rank_data(sushi_rankings[i, ])
    )
  }
  expect_equal(
    mean(mod$alpha_samples),
    mean(mod_bmm$alpha_samples[-(1:500)]),
    tolerance = .01
  )
  expect_equal(
    sd(mod$alpha_samples),
    sd(mod_bmm$alpha_samples[-(1:500)]),
    tolerance = .04
  )

  dat <- sushi_rankings[sample(nrow(sushi_rankings), n), ]
  dat[dat > 6] <- NA
  mod_bmm <- compute_mallows(data = setup_rank_data(dat))

  mod <- sample_prior(1000, ncol(dat))
  for (i in seq_len(n)) {
    mod <- update_mallows(
      model = mod,
      new_data = setup_rank_data(dat[i, , drop = FALSE])
    )
  }
  expect_equal(
    mean(mod$alpha_samples),
    mean(mod_bmm$alpha_samples[-(1:500)]),
    tolerance = .02
  )
  expect_equal(
    sd(mod$alpha_samples),
    sd(mod_bmm$alpha_samples[-(1:500)]),
    tolerance = .1
  )
})

test_that("update_mallows works with swap prosal for rho", {
  triple_potato <- rbind(potato_visual, potato_visual, potato_visual)
  rownames(triple_potato) <- seq_len(nrow(triple_potato))
  user_ids <- rownames(triple_potato)

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(triple_potato),
    compute_options = set_compute_options(
      nmc = 10000, burnin = 1000, rho_proposal = "swap")
  )

  mod_init <- compute_mallows(
    data = setup_rank_data(triple_potato[1:4, , drop = FALSE]),
    compute_options = set_compute_options(
      nmc = 20000, burnin = 1000, rho_proposal = "swap")
  )

  mod_smc <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = triple_potato[5:20, ])
  )

  mod_smc_next <- update_mallows(
    model = mod_smc,
    new_data = setup_rank_data(rankings = triple_potato[21:36, ])
  )

  expect_equal(
    mean(mod_smc_next$alpha$value),
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    tolerance = 0.01
  )

  expect_equal(
    sd(mod_smc_next$alpha$value),
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    tolerance = 0.1
  )
})
