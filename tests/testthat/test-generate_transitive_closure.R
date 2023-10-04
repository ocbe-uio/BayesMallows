test_that("generate_transitive_closure works", {
  pair_comp <- data.frame(
    assessor = c(1, 1, 1, 2, 2),
    bottom_item = c(2, 5, 3, 5, 3),
    top_item = c(1, 1, 5, 3, 4)
  )

  # We then generate the transitive closure of these preferences:
  res1 <- generate_transitive_closure(pair_comp)

  cl <- parallel::makeCluster(2)
  res2 <- generate_transitive_closure(pair_comp)
  parallel::stopCluster(cl)

  expect_equal(
    res1,
    structure(list(assessor = c(1, 1, 1, 1, 2, 2, 2), bottom_item = c(
      2,
      3, 5, 3, 5, 3, 5
    ), top_item = c(1, 1, 1, 5, 3, 4, 4)), row.names = c(
      NA,
      -7L
    ), class = c("BayesMallowsTC", "data.frame"))
  )
  expect_equal(res1, res2)

  pair_comp <- rbind(pair_comp, data.frame(assessor = 1, bottom_item = 1, top_item = 3))

  # This causes an error message and prints out the problematic rankings:
  expect_error(
    capture.output(generate_transitive_closure(pair_comp)),
    "Cannot compute transitive closure. Please run compute_mallows with error_model='bernoulli'."
  )
})
