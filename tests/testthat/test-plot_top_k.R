test_that("predict_top_k works", {
  set.seed(1)
  model_fit <- compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(
      nmc = 100, burnin = 0, save_aug = TRUE
    )
  )

  ptk <- predict_top_k(model_fit)
  expect_equal(ptk$prob[[900]], .06)
  expect_equal(dim(ptk), c(900, 3))

  p <- plot_top_k(model_fit)
  expect_equal(dim(p$data), c(900, 3))
  expect_equal(p$labels$fill, "Prob.")
  expect_equal(p$labels$y, "Item")
  expect_equal(p$labels$x, "Assessor")
  expect_s3_class(p, "ggplot")

  model_fit <- compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 10, burnin = 0)
  )

  expect_error(predict_top_k(model_fit), "Please set save_aug = TRUE")
  expect_error(plot_top_k(model_fit), "Please set save_aug = TRUE")

  model_fit <- compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 10)
  )

  expect_error(predict_top_k(model_fit), "Please specify the burnin.")
})
