test_that("Mixtures deals properly with thinning", {
  bmm <- compute_mallows(
    data = setup_rank_data(potato_visual),
    model_options = set_model_options(n_cluster = 3),
    compute_options = set_compute_options(
      nmc = 100, burnin = 50,
      clus_thinning = 20, rho_thinning = 20, include_wcd = TRUE)
  )

  expect_s3_class(assign_cluster(bmm), "data.frame")
})
