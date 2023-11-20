test_that("generate_transitive_closure works", {
  pair_comp <- data.frame(
    assessor = c(1, 1, 1, 2, 2),
    bottom_item = c(2, 5, 3, 5, 3),
    top_item = c(1, 1, 5, 3, 4)
  )

  dat <- setup_rank_data(preferences = pair_comp)

  cl <- parallel::makeCluster(2)
  dat2 <- setup_rank_data(preferences = pair_comp, cl = cl)
  parallel::stopCluster(cl)

  expect_equal(
    dat$preferences,
    structure(list(assessor = c(1, 1, 1, 1, 2, 2, 2), bottom_item = c(2, 3, 5, 3, 5, 3, 5), top_item = c(
      1, 1, 1,
      5, 3, 4, 4
    )), row.names = c(NA, -7L), class = c(
      "BayesMallowsTransitiveClosure",
      "data.frame"
    ))
  )
  expect_equal(dat$preferences, dat2$preferences)

  pair_comp <- rbind(pair_comp, data.frame(assessor = 1, bottom_item = 1, top_item = 3))

  # This causes an error message and prints out the problematic rankings:
  expect_message(
    setup_rank_data(preferences = pair_comp),
    "Preferences are intransitive."
  )
})
