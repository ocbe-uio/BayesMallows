test_that("predict_top_k works", {
  set.seed(1)
  model_fit <- compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(
      nmc = 1000, burnin = 0, save_aug = TRUE))

  tail(predict_top_k(model_fit))
})
