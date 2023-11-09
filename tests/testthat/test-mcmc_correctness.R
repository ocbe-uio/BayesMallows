test_that("compute_mallows is correct for complete data", {

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  expect_equal(
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    10.85, tolerance = .05)

  expect_equal(
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    0.77, tolerance = .1)

})


test_that("compute_mallows is correct for pairwise preferences", {

  set.seed(123)

  mod_bmm <- compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  expect_equal(
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    4.825, tolerance = .05)

  expect_equal(
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    0.286, tolerance = .1)

})


test_that("compute_mallows is correct for top-k ranks", {

  set.seed(123)
  dat <- potato_visual
  dat[dat > 10] <- NA

  mod_bmm <- compute_mallows(
    data = setup_rank_data(dat),
    compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

  expect_equal(
    mean(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    10.29, tolerance = .01)

  expect_equal(
    sd(mod_bmm$alpha$value[mod_bmm$alpha$iteration > 1000]),
    1.506, tolerance = .1)

})
