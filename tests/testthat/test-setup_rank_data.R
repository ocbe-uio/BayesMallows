test_that("setup_rank_data works with rankings", {
  expect_error(
    setup_rank_data(),
    "Either rankings or preferences")

  rr <- potato_visual
  rr[1, 1] <- NA
  expect_error(
    setup_rank_data(rr, na_action = "fail"),
   "rankings matrix contains NA values")

  expect_snapshot(
    dat <- setup_rank_data(rr, na_action = "omit"))

  expect_equal(dim(dat$rankings), c(11, 20))

  dat <- setup_rank_data(potato_visual,
                         observation_frequency = 1:12)
  dat2 <- setup_rank_data(potato_visual)
  expect_equal(dat$rankings, dat2$rankings)

  expect_error(
    setup_rank_data(
      rankings = potato_visual,
      observation_frequency = 1:19),
    "observation_frequency must be of same length as the number of rows in rankings"
  )

  rr <- matrix(c(1, 1, 2, 1), ncol = 2)
  expect_error(
    setup_rank_data(rr),
    "invalid permutations provided in rankings matrix")

  dat <- setup_rank_data(rr, validate_rankings = FALSE)
  expect_equal(dat$rankings, rr)
})
