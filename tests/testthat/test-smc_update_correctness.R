test_that("update_mallows is correct for new rankings", {

  triple_potato <- rbind(potato_visual, potato_visual, potato_visual)

  set.seed(123)

  mod_bmm <- compute_mallows(rankings = triple_potato, nmc = 10000)
  mod_bmm$burnin <- 1000

  mod_init <- compute_mallows(rankings = triple_potato[1:4, , drop = FALSE],
                              nmc = 10000)
  mod_init$burnin <- 1000
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
               tolerance = 0.02)

  expect_equal(mean(mod_smc_next$alpha$value), 10.76, tolerance = .1)

  expect_equal(sd(mod_smc_next$alpha$value),
               sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
               tolerance = 0.02)

  expect_equal(sd(mod_smc_next$alpha$value), .43, tolerance = .01)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  # How many items are in disagreement
  expect_lte(
    rank_distance(
      matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
      as.numeric(as.factor(smc_consensus$item)),
      metric = "ulam"
    ),
    1
  )

  set.seed(123)

  mod_bmm <- compute_mallows(rankings = sushi_rankings, nmc = 1000)
  mod_bmm$burnin <- 200

  mod_init <- compute_mallows(rankings = sushi_rankings[1:100, ],
                              nmc = 10000)
  mod_init$burnin <- 1000
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
               tolerance = 0.02)

  expect_equal(mean(mod_smc_next$alpha$value), 1.7, tolerance = 0.02)

  expect_equal(sd(mod_smc_next$alpha$value),
               sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 200]),
               tolerance = 0.01)

  expect_equal(sd(mod_smc_next$alpha$value), 0.01, tolerance = 0.005)

  bmm_consensus <- compute_consensus(mod_bmm)
  smc_consensus <- compute_consensus(mod_smc_next)

  # How many items are in disagreement
  expect_lte(
    rank_distance(
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

  bmm_mod <- compute_mallows(rankings = dat, nmc = 10000)
  bmm_mod$burnin <- 1000

  mod_init <- compute_mallows(
    rankings = dat[1:4, , drop = FALSE], nmc = 10000, save_aug = TRUE)
  mod_init$burnin <- 1000

  for(aug in c("uniform", "pseudo")) {
    mod_smc <- update_mallows(
      model = mod_init,
      new_rankings = dat[5:20, ],
      n_particles = 3000,
      mcmc_steps = 5,
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
      tolerance = ifelse(aug == "uniform", .2, .05)
    )

    bmm_consensus <- compute_consensus(bmm_mod)
    smc_consensus <- compute_consensus(mod_smc_next)

    expect_lte(
      rank_distance(
        matrix(as.numeric(as.factor(bmm_consensus$item)), nrow = 1),
        as.numeric(as.factor(smc_consensus$item)),
        metric = "ulam"
      ),
      ifelse(aug == "uniform", 4, 1)
    )
  }

})
