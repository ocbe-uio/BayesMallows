beach_small <- subset(
  beach_preferences,
  bottom_item %in% 1:5 & top_item %in% 1:5
)
dat <- setup_rank_data(preferences = beach_small)

test_that("plot_top_k and predict_top_k fail when they should", {
  bmm <- compute_mallows(
    data = dat,
    compute_options = set_compute_options(nmc = 2, save_aug = TRUE)
  )
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

  bmm <- compute_mallows(
    data = dat,
    compute_options = set_compute_options(nmc = 2, save_aug = FALSE)
  )
  # Expecting error because save_aug = FALSE
  expect_error(plot_top_k(bmm, burnin = 50))
  expect_error(predict_top_k(bmm, burnin = 50))

  # Test whether predict_top_k returns correct numbers
  set.seed(1)
  bmm <- compute_mallows(
    data = dat,
    compute_options = set_compute_options(nmc = 20, save_aug = TRUE)
  )

  pred <- predict_top_k(bmm, burnin = 4, k = 3)

  expect_equal(
    pred[order(pred$assessor, pred$item), ][13:16, ],
    structure(
      list(
        assessor = c(3, 3, 3, 4), item = c("Item 3", "Item 4", "Item 5", "Item 1"),
        prob = c(1, 0.3125, 0, 0.5625)
      ),
      row.names = c(119L, 177L, 235L, 4L), class = "data.frame"
    )
  )

  pred <- predict_top_k(bmm, burnin = 4, k = 5)
  expect_equal(
    head(pred),
    structure(list(
      assessor = c(1, 1, 1, 1, 1, 2),
      item = c("Item 1", "Item 2", "Item 3", "Item 4", "Item 5", "Item 1"),
      prob = c(1, 1, 1, 1, 1, 1)
    ), row.names = c(1L, 59L, 117L, 175L, 233L, 2L), class = "data.frame")
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
  set.seed(1)
  bmm1 <- compute_mallows(
    setup_rank_data(dat),
    compute_options = set_compute_options(nmc = 500, save_aug = TRUE, aug_thinning = 1L)
  )
  set.seed(1)
  bmm3 <- compute_mallows(
    setup_rank_data(dat),
    compute_options = set_compute_options(nmc = 500, save_aug = TRUE, aug_thinning = 3L)
  )

  pred1 <- predict_top_k(bmm1, burnin = 200, k = k)
  pred3 <- predict_top_k(bmm3, burnin = 200, k = k)

  expect_equal(pred1, pred3)
})

dat <- setup_rank_data(preferences = beach_preferences)
test_that("plot_top_k works", {
  bmm <- compute_mallows(
    dat,
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE)
  )
  expect_s3_class(plot_top_k(bmm, k = 4, burnin = 5), "ggplot")

  bmm <- compute_mallows(
    dat,
    model = set_model_options(n_clusters = 2),
    compute_options = set_compute_options(nmc = 10, save_aug = TRUE)
  )
  expect_s3_class(plot_top_k(bmm, k = 4, burnin = 5, rel_widths = c(.5, 1)), "ggplot")
})
