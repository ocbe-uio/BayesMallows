test_that("setup_rank_data works with rankings", {
  expect_error(
    setup_rank_data(),
    "Either rankings or preferences"
  )

  rr <- potato_visual
  rr[1, 1] <- NA
  expect_error(
    setup_rank_data(rr, na_action = "fail"),
    "rankings matrix contains NA values"
  )

  expect_snapshot(
    dat <- setup_rank_data(rr, na_action = "omit")
  )

  expect_equal(dim(dat$rankings), c(11, 20))

  dat <- setup_rank_data(potato_visual,
    observation_frequency = 1:12
  )
  dat2 <- setup_rank_data(potato_visual)
  expect_equal(dat$rankings, dat2$rankings)

  expect_error(
    setup_rank_data(
      rankings = potato_visual,
      observation_frequency = 1:19
    ),
    "observation_frequency must be of same length as the number of rows in rankings"
  )

  rr <- matrix(c(1, 1, 2, 1), ncol = 2)
  expect_error(
    setup_rank_data(rr),
    "invalid permutations provided in rankings matrix"
  )

  dat <- setup_rank_data(rr, validate_rankings = FALSE)
  expect_equal(dat$rankings, rr)
})

test_that("setup_rank_data works for preferences", {
  expect_error(
    setup_rank_data(preferences = beach_preferences, random = TRUE),
    "Number of items exceeds the limit"
  )

  rr <- matrix(rep(1:3, 2), byrow = TRUE, ncol = 3)
  pp <- data.frame(assessor = 1:2, bottom_item = 2, top_item = 1)
  dat <- setup_rank_data(rr, pp)
  expect_equal(dat$rankings, rr)
  expect_equal(as.data.frame(dat$preferences), pp)

  rr <- matrix(rep(1:2, 2), byrow = TRUE, ncol = 2)
  pp <- data.frame(assessor = 1:2, bottom_item = 3, top_item = 1)
  dat <- setup_rank_data(rankings = rr, preferences = pp)

  cl <- parallel::makeCluster(2)
  dat2 <- setup_rank_data(rankings = rr, preferences = pp, cl = cl)
  parallel::stopCluster(cl)

  expect_equal(dat2$rankings, dat$rankings)
  expect_equal(dat2$preferences, dat$preferences)

  set.seed(1)
  prefdat <- subset(beach_preferences, assessor <= 3)
  dat <- setup_rank_data(
    preferences = prefdat, shuffle_unranked = TRUE
  )
  expect_equal(dim(dat$rankings), c(3, 15))
  expect_equal(sum(is.na(dat$rankings)), 0)
  dat2 <- setup_rank_data(
    preferences = prefdat, shuffle_unranked = TRUE
  )
  expect_false(
    all(dat2$rankings == dat$rankings)
  )
  expect_equal(dat2$preferences, dat$preferences)

  prefdat <- subset(
    beach_preferences,
    bottom_item <= 3 & top_item <= 3
  )

  dat1 <- setup_rank_data(preferences = prefdat, random = TRUE)
  expect_error(
    setup_rank_data(preferences = prefdat, random = TRUE, random_limit = 1),
    "Number of items exceeds"
  )
  dat2 <- setup_rank_data(
    preferences = prefdat, random = TRUE,
    shuffle_unranked = TRUE
  )

  expect_false(all(dat1$rankings == dat2$rankings))
  expect_equal(dim(dat1$rankings), dim(dat2$rankings))

  prefdat$assessor <- as.character(prefdat$assessor)
  expect_error(
    setup_rank_data(preferences = prefdat),
    "assessor column in preferences must be numeric"
  )
})
