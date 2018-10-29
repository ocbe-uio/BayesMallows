context("Testing plot_top_k and predict_top_k")

test_that("plot_top_k and predict_top_k fail when they should", {
  beach_tc <- generate_transitive_closure(beach_preferences)
  beach_init_rank <- generate_initial_ranking(beach_tc)
  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                         nmc = 100, save_aug = TRUE)
  # Expecting error because burnin is not set
  expect_error(plot_top_k(bmm))
  expect_error(predict_top_k(bmm))

  # Expecting error because burnin > nmc
  expect_error(plot_top_k(bmm, burnin = 200))
  expect_error(plot_top_k(bmm, burnin = 100))
  bmm$burnin <- 100
  expect_error(plot_top_k(bmm))
  bmm$burnin <- "text"
  expect_error(plot_top_k(bmm))

  expect_error(predict_top_k(bmm, burnin = 200))
  expect_error(predict_top_k(bmm, burnin = 100))
  bmm$burnin <- 100
  expect_error(predict_top_k(bmm))
  bmm$burnin <- "text"
  expect_error(predict_top_k(bmm))

  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                         nmc = 100, save_aug = FALSE)
  # Expecting error because save_aug = FALSE
  expect_error(plot_top_k(bmm, burnin = 50))
  expect_error(predict_top_k(bmm, burnin = 50))
})
