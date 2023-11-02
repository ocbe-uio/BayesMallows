
test_that("skip_extended works", {
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "true")
  expect_condition(skip_extended(), NA, class = "skip")
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "")
  expect_condition(skip_extended(), class = "skip")
})


test_that("update_mallows is correct for new rankings", {

  triple_potato <- rbind(potato_visual, potato_visual, potato_visual)

  set.seed(123)

  # Metropolis-Hastings
  mod_bmm <- compute_mallows(rankings = triple_potato, nmc = 10000)
  mod_bmm$burnin <- 1000

  # Sequentially, using update function
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

  # Metropolis-Hastings
  mod_bmm <- compute_mallows(rankings = sushi_rankings, nmc = 1000)
  mod_bmm$burnin <- 200

  # Sequentially, using update function
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
  skip_extended()

  # Small dummy data ----
  set.seed(123)

  rankings <- matrix(rep(c(
    1, 2, 3,
    1, 3, 2,
    1, 2, 3,
    1, 2, 3,
    2, 1, 3
  ), times = 10), ncol = 3, byrow = TRUE)

  rankings[sample(seq_along(rankings), 10)] <- NA

  bmm_mod <- compute_mallows(rankings = rankings, nmc = 10000)
  bmm_mod$burnin <- 1000

  inds <- rep(1:10, each = 5)
  smc_init <- compute_mallows(rankings = rankings[inds == 1, ],
                              nmc = 10000)
  smc_init$burnin <- 5000

  smc_update <- smc_init
  for (i in 2:10) {
    smc_update <- update_mallows(
      model = smc_update,
      new_rankings = rankings[inds == i, ],
      type = "partial",
      n_particles = 10000,
      verbose = TRUE
    )
  }

  expect_equal(mean(smc_update$alpha$value), 2.51, tolerance = .01)
  expect_equal(mean(bmm_mod$alpha$value[bmm_mod$alpha$iteration > 100]), 2.48,
    tolerance = .01
  )

  expect_equal(consensus_unwrapper(smc_update), 1:3)
  expect_equal(consensus_unwrapper(bmm_mod), 1:3)

  # Sushi data with missingness ----
  dat <- sushi_rankings
  dat[dat > 8] <- NA
  bmm_sushi <- compute_mallows(rankings = dat)
  bmm_sushi$burnin <- 300

  smc_sushi_init <- smc_mallows_new_users(
    rankings = dat[1:1000, ],
    type = "partial",
    n_particles = 1000,
    timesteps = 1,
    num_new_obs = 1000,
    mcmc_steps = 5
  )

  smc_sushi_update <- smc_sushi_init
  inds <- rep(1:5, each = 1000)
  for (i in 2:5) {
    smc_sushi_update <- smc_mallows_update(
      model = smc_sushi_update, rankings = dat[inds == i, ]
    )
  }

  expect_equal(mean(bmm_sushi$alpha$value[bmm_sushi$alpha$iteration > 300]),
    1.69,
    tolerance = .01
  )
  expect_equal(mean(smc_sushi_update$alpha_samples[, 1]),
    1.73,
    tolerance = .01
  )

  smc_consensus <- compute_consensus(smc_sushi_update)
  smc_consensus$item <- unlist(
    regmatches(
      smc_consensus$item,
      gregexpr("[0-9]+", smc_consensus$item)
    )
  )
  smc_consensus$item <- colnames(sushi_rankings)[as.integer(smc_consensus$item)]
  bmm_consensus <- compute_consensus(bmm_sushi)

  expect_equal(
    bmm_consensus$ranking,
    smc_consensus$ranking
  )
})
