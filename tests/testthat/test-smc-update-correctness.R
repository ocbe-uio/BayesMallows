consensus_unwrapper <- function(x) {
  x <- compute_consensus(x)
  as.numeric(unlist(regmatches(x$item, gregexpr("[0-9]+", x$item))))
}


test_that("skip_extended works", {
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "true")
  expect_condition(skip_extended(), NA, class = "skip")
  Sys.setenv(BAYESMALLOWS_EXTENDED_TESTS = "")
  expect_condition(skip_extended(), class = "skip")
})


test_that("update_mallows is correct for new rankings", {

  set.seed(123)

  # Metropolis-Hastings
  mod_bmm <- compute_mallows(rankings = potato_visual, nmc = 10000)
  mod_bmm$burnin <- 1000

  # Sequentially, using update function
  mod_init <- compute_mallows(rankings = potato_visual[1:6, ], nmc = 10000)
  mod_init$burnin <- 1000
  mod_smc <- update_mallows(
    model = mod_init,
    new_rankings = potato_visual[7:12, ],
    n_particles = 500
    )


  # Poterior mean of alpha should be the same in both SMC methods, and close to BMM
  expect_equal(mean(mod_smc$alpha$value),
               mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
               tolerance = 0.1)

  # Is there any disagreement between the methods about the ranking of the items?
  bmm_consensus <- consensus_unwrapper(mod_bmm)
  smc_consensus <- consensus_unwrapper(mod_smc)

  # How many items are in disagreement
  expect_equal(
    rank_distance(matrix(bmm_consensus, nrow = 1),
      smc_consensus,
      metric = "ulam"
    ),
    1
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
