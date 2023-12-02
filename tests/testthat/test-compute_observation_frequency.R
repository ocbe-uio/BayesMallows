test_that("compute_observation_frequency works", {
  expect_error(
    compute_observation_frequency(rankings = beach_preferences)
  )

  mat <- matrix(c(rep(1:3, 5), rep(3:1, 10), rep(c(1, 3, 2), 2)),
    byrow = TRUE, ncol = 3
  )

  freq <- compute_observation_frequency(rankings = mat)[, 4]
  expect_equal(freq, c(5, 2, 10))

  mat[mat < 3] <- NA
  freq <- compute_observation_frequency(mat)[, 4]
  expect_equal(freq, c(5, 2, 10))
})
