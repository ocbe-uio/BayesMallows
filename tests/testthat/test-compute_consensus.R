test_that("compute_consensus fails properly", {
  mod <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10)
  )
  expect_error(compute_consensus(mod), "Please specify the burnin")
  mod$burnin <- 11
  expect_error(compute_consensus(mod), "burnin < model_fit")
  mod$burnin <- 2
  expect_error(
    compute_consensus(mod, parameter = "Rtilde"),
    "For augmented ranks, please refit model"
  )

  dat <- potato_visual
  dat[c(1, 13, 14, 23)] <- NA
  mod <- compute_mallows(
    setup_rank_data(dat),
    compute_options = set_compute_options(nmc = 10, burnin = 2)
  )

  expect_error(
    compute_consensus(mod, parameter = "Rtilde"),
    "For augmented ranks, please refit"
  )
})

test_that("compute_consensus.BayesMallows works", {
  set.seed(1)
  mod <- compute_mallows(
    setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 200, burnin = 100)
  )
  c1 <- compute_consensus(mod)
  expect_gt(which(c1$item == "P16") - which(c1$item == "P12"), 0)
  c2 <- compute_consensus(mod, type = "MAP")
  expect_equal(length(unique(c2$probability)), 1)

  mod <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 100, burnin = 50, save_aug = TRUE)
  )
  a1 <- compute_consensus(mod, parameter = "Rtilde", assessors = 2)
  expect_equal(unique(a1$assessor), 2)
  expect_equal(dim(a1), c(15, 4))
  a2 <- compute_consensus(
    mod,
    parameter = "Rtilde", type = "MAP", assessors = 3
  )
  expect_equal(unique(a2$assessor), 3)
  expect_equal(length(unique(a2$probability)), 1)
})

test_that("compute_consensus.SMCMallows works", {
  set.seed(1)
  data_first_batch <- potato_visual[1:4, ]
  mod_init <- compute_mallows(
    data = setup_rank_data(data_first_batch),
    compute_options = set_compute_options(nmc = 200, burnin = 50)
  )

  data_second_batch <- potato_visual[5:8, ]
  mod_second <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = data_second_batch),
    smc_options = set_smc_options(n_particles = 30)
  )

  expect_equal(dim(compute_consensus(mod_second)), c(20, 4))

  data_third_batch <- potato_visual[9:12, ]
  mod_final <- update_mallows(
    model = mod_second,
    new_data = setup_rank_data(rankings = data_third_batch)
  )

  expect_error(
    compute_consensus(mod_final, parameter = "Rtilde"),
    "'arg' should be"
  )

  a1 <- compute_consensus(mod_final, type = "MAP")
  expect_equal(length(unique(a1$probability)), 1)
  a2 <- compute_consensus(mod_final, type = "CP")
  expect_equal(dim(a2), c(20, 4))

  mod <- sample_prior(
    1000, ncol(sushi_rankings),
    priors = set_priors(gamma = 2, lambda = .1))
  for (i in seq_len(30)) {
    mod <- update_mallows(
      model = mod,
      new_data = setup_rank_data(sushi_rankings[i, , drop = FALSE])
    )
  }
  expect_equal(
    compute_consensus(mod)$item,
    c("fatty tuna", "shrimp", "squid", "salmon roe", "tuna", "sea urchin",
      "tuna roll", "sea eel", "egg", "cucumber roll")
  )
})
