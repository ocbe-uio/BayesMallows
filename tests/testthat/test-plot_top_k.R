context("Testing plot_top_k and predict_top_k")
library(dplyr)
beach_small <- beach_preferences %>%
  filter(bottom_item %in% 1:5, top_item %in% 1:5)
beach_tc <- generate_transitive_closure(beach_small)
beach_init_rank <- generate_initial_ranking(beach_tc)

test_that("plot_top_k and predict_top_k fail when they should", {

  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                         nmc = 2, save_aug = TRUE)
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
                         nmc = 2, save_aug = FALSE)
  # Expecting error because save_aug = FALSE
  expect_error(plot_top_k(bmm, burnin = 50))
  expect_error(predict_top_k(bmm, burnin = 50))

  # Test whether predict_top_k returns correct numbers
  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                         nmc = 20, save_aug = TRUE, seed = 1L)

  pred <- predict_top_k(bmm, burnin = 4, k = 3)
  expect_equal(
    pred[13:16, ],
    structure(list(assessor = c(3, 3, 3, 4), item = c("Item 3", "Item 4",
                                                      "Item 5", "Item 1"), prob = c(1, 0.3125, 0, 0.5625)), row.names = c("3.3",
                                                                                                                          "3.4", "3.5", "4.1"), class = "data.frame")
    )

  pred <- predict_top_k(bmm, burnin = 4, k = 5)
  expect_equal(
    head(pred),
    structure(list(assessor = c(1, 1, 1, 1, 1, 2), item = c("Item 1",
                                                            "Item 2", "Item 3", "Item 4", "Item 5", "Item 1"), prob = c(1,
                                                                                                                        1, 1, 1, 1, 1)), row.names = c("1.1", "1.2", "1.3", "1.4", "1.5",
                                                                                                                                                       "2.1"), class = "data.frame")
    )
})


test_that("predict_top_k works with augmentation thinning", {
  # Create top-k dataset
  k <- 3
  dat <- potato_visual
  dat[dat > 3] <- NA
  dat <- dat[, apply(dat, 2, function(x) !all(is.na(x)))]

  # Run the model with different levels of thinning
  # This model is deterministic, because each subject has either ranked or not ranked the item
  bmm1 <- compute_mallows(rankings = dat, nmc = 500, save_aug = TRUE, aug_thinning = 1L,
                          seed = 1L)
  bmm3 <- compute_mallows(rankings = dat, nmc = 500, save_aug = TRUE, aug_thinning = 3L,
                          seed = 1L)

  pred1 <- predict_top_k(bmm1, burnin = 200, k = k)
  pred3 <- predict_top_k(bmm3, burnin = 200, k = k)

  expect_equal(pred1, pred3)
})


test_that("plot_top_k works", {
  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                         nmc = 10, save_aug = TRUE)
  expect_s3_class(plot_top_k(bmm, k = 4, burnin = 5), "ggplot")

  bmm <- compute_mallows(rankings = beach_init_rank, preferences = beach_tc,
                               nmc = 10, save_aug = TRUE, n_clusters = 2)
  expect_s3_class(plot_top_k(bmm, k = 4, burnin = 5, rel_widths = c(.5, 1)), "ggplot")

})
