context("Testing compute_consensus")

b <- compute_mallows(preferences = beach_preferences, nmc = 500, seed = 123L)
b$burnin <- 200
cp <- compute_consensus(b)
map <- compute_consensus(b, type = "MAP")

test_that("compute_consensus returns correct object", {
  expect_true(inherits(cp, "data.frame"))
  expect_true(inherits(map, "data.frame"))

})

test_that("compute_consensus computes correctly for rho", {
  expect_equal(cp, structure(list(ranking = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                              12, 13, 14, 15), item = c("Item 9", "Item 6", "Item 3", "Item 11",
                                                                        "Item 15", "Item 10", "Item 1", "Item 5", "Item 7", "Item 13",
                                                                        "Item 8", "Item 4", "Item 12", "Item 14", "Item 2"), cumprob = c(0.78,
                                                                                                                                         1, 1, 0.976666666666667, 0.87, 0.983333333333333, 1, 0.703333333333333,
                                                                                                                                         1, 1, 0.926666666666667, 0.926666666666667, 1, 0.913333333333333,
                                                                                                                                         1)), row.names = c(NA, -15L), class = c("tbl_df", "tbl", "data.frame"
                                                                                                                                         )))

  expect_equal(map, structure(list(probability = c(0.34, 0.34, 0.34, 0.34, 0.34,
                                                   0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34),
                                   item = c("Item 9", "Item 6", "Item 3", "Item 11", "Item 15",
                                            "Item 10", "Item 1", "Item 5", "Item 7", "Item 13", "Item 8",
                                            "Item 4", "Item 12", "Item 14", "Item 2"), map_ranking = c(1,
                                                                                                       2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)), row.names = c(NA,
                                                                                                                                                                       -15L), class = c("tbl_df", "tbl", "data.frame")))


})


b2 <- compute_mallows(preferences = beach_preferences, nmc = 500, save_aug = TRUE, seed = 123L)
b3 <- compute_mallows(preferences = beach_preferences, nmc = 500, save_aug = TRUE, aug_thinning = 3, seed = 123L)

test_that("compute_consensus fails when it should", {
  # Burnin not set
  expect_error(compute_consensus(b2, type = "CP"))
  expect_error(compute_consensus(b2, type = "MAP"))

  # Augmented data are missing
  expect_error(compute_consensus(b, burnin = 200, type = "CP", parameter = "Rtilde"))
  expect_error(compute_consensus(b, burnin = 200, type = "MAP", parameter = "Rtilde"))

  # Incorrect assessor
  expect_error(compute_consensus(b2, burnin = 200, type = "CP", parameter = "Rtilde", assessors = 0))
  expect_error(compute_consensus(b2, burnin = 200, type = "MAP", parameter = "Rtilde", assessors = 0))

})

test_that("compute_consensus computes augmented ranks correctly", {
  res <- compute_consensus(b2, type = "CP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(ranking = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                               12, 13, 14, 15), item = c("Item 6", "Item 11", "Item 3", "Item 9",
                                                                         "Item 10", "Item 15", "Item 1", "Item 12", "Item 7", "Item 13",
                                                                         "Item 5", "Item 4", "Item 8", "Item 14", "Item 2"), cumprob = c(0.52,
                                                                                                                                         0.903333333333333, 0.696666666666667, 0.46, 0.683333333333333,
                                                                                                                                         0.623333333333333, 0.666666666666667, 0.76, 0.67, 0.583333333333333,
                                                                                                                                         0.433333333333333, 0.53, 0.76, 0.626666666666667, 1)), row.names = c(NA,
                                                                                                                                                                                                              -15L), class = c("tbl_df", "tbl", "data.frame")))


  res <- compute_consensus(b3, type = "CP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(ranking = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                               12, 13, 14, 15), item = c("Item 6", "Item 11", "Item 3", "Item 9",
                                                                         "Item 10", "Item 15", "Item 1", "Item 12", "Item 7", "Item 13",
                                                                         "Item 5", "Item 4", "Item 8", "Item 14", "Item 2"), cumprob = c(0.52,
                                                                                                                                         0.91, 0.7, 0.46, 0.69, 0.63, 0.67, 0.75, 0.67, 0.58, 0.44, 0.54,
                                                                                                                                         0.74, 0.62, 1)), row.names = c(NA, -15L), class = c("tbl_df",
                                                                                                                                                                                             "tbl", "data.frame")))



  res <- compute_consensus(b2, type = "CP", burnin = 200, parameter = "Rtilde", assessors = c(3L, 5L))
  expect_equal(res, structure(list(assessor = c("3", "3", "3", "3", "3", "3", "3",
                                                "3", "3", "3", "3", "3", "3", "3", "3", "5", "5", "5", "5", "5",
                                                "5", "5", "5", "5", "5", "5", "5", "5", "5", "5"), ranking = c(1,
                                                                                                               2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5,
                                                                                                               6, 7, 8, 9, 10, 11, 12, 13, 14, 15), item = c("Item 9", "Item 6",
                                                                                                                                                             "Item 15", "Item 1", "Item 3", "Item 11", "Item 10", "Item 7",
                                                                                                                                                             "Item 13", "Item 14", "Item 8", "Item 5", "Item 2", "Item 12",
                                                                                                                                                             "Item 4", "Item 9", "Item 6", "Item 10", "Item 11", "Item 5",
                                                                                                                                                             "Item 4", "Item 3", "Item 14", "Item 15", "Item 1", "Item 7",
                                                                                                                                                             "Item 8", "Item 2", "Item 12", "Item 13"), cumprob = c(1, 0.916666666666667,
                                                                                                                                                                                                                    0.733333333333333, 0.363333333333333, 0.506666666666667, 0.57,
                                                                                                                                                                                                                    0.65, 0.513333333333333, 0.69, 0.6, 0.586666666666667, 0.77,
                                                                                                                                                                                                                    0.64, 0.67, 1, 0.39, 0.783333333333333, 0.786666666666667, 0.703333333333333,
                                                                                                                                                                                                                    0.76, 0.73, 0.57, 0.326666666666667, 0.486666666666667, 0.68,
                                                                                                                                                                                                                    0.78, 0.853333333333333, 0.596666666666667, 0.833333333333333,
                                                                                                                                                                                                                    1)), row.names = c(NA, -30L), class = c("tbl_df", "tbl", "data.frame"
                                                                                                                                                                                                                    )))


  res <- compute_consensus(b3, type = "CP", burnin = 200, parameter = "Rtilde", assessors = c(3L, 5L))
  expect_equal(res, structure(list(assessor = c("3", "3", "3", "3", "3", "3", "3",
                                                "3", "3", "3", "3", "3", "3", "3", "3", "5", "5", "5", "5", "5",
                                                "5", "5", "5", "5", "5", "5", "5", "5", "5", "5"), ranking = c(1,
                                                                                                               2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5,
                                                                                                               6, 7, 8, 9, 10, 11, 12, 13, 14, 15), item = c("Item 9", "Item 6",
                                                                                                                                                             "Item 15", "Item 1", "Item 3", "Item 11", "Item 10", "Item 7",
                                                                                                                                                             "Item 13", "Item 14", "Item 8", "Item 5", "Item 2", "Item 12",
                                                                                                                                                             "Item 4", "Item 9", "Item 6", "Item 10", "Item 11", "Item 5",
                                                                                                                                                             "Item 4", "Item 3", "Item 14", "Item 15", "Item 1", "Item 7",
                                                                                                                                                             "Item 8", "Item 2", "Item 12", "Item 13"), cumprob = c(1, 0.92,
                                                                                                                                                                                                                    0.75, 0.36, 0.52, 0.56, 0.64, 0.52, 0.7, 0.61, 0.57, 0.78, 0.64,
                                                                                                                                                                                                                    0.68, 1, 0.38, 0.79, 0.78, 0.71, 0.75, 0.73, 0.57, 0.33, 0.48,
                                                                                                                                                                                                                    0.66, 0.78, 0.85, 0.6, 0.84, 1)), row.names = c(NA, -30L), class = c("tbl_df",
                                                                                                                                                                                                                                                                                         "tbl", "data.frame")))



  res <- compute_consensus(b2, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(probability = c(0.0433333333333333, 0.0433333333333333,
                                                   0.0433333333333333, 0.0433333333333333, 0.0433333333333333, 0.0433333333333333,
                                                   0.0433333333333333, 0.0433333333333333, 0.0433333333333333, 0.0433333333333333,
                                                   0.0433333333333333, 0.0433333333333333, 0.0433333333333333, 0.0433333333333333,
                                                   0.0433333333333333), item = c("Item 6", "Item 11", "Item 3",
                                                                                 "Item 9", "Item 10", "Item 15", "Item 1", "Item 12", "Item 7",
                                                                                 "Item 13", "Item 4", "Item 5", "Item 8", "Item 2", "Item 14"),
                                   map_ranking = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                                   14, 15)), row.names = c(NA, -15L), class = c("tbl_df", "tbl",
                                                                                                "data.frame")))


  res <- compute_consensus(b3, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(probability = c(0.05, 0.05, 0.05, 0.05, 0.05,
                                                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
                                   item = c("Item 6", "Item 11", "Item 3", "Item 9", "Item 10",
                                            "Item 15", "Item 1", "Item 12", "Item 7", "Item 13", "Item 4",
                                            "Item 5", "Item 8", "Item 2", "Item 14"), map_ranking = c(1,
                                                                                                      2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)), row.names = c(NA,
                                                                                                                                                                      -15L), class = c("tbl_df", "tbl", "data.frame")))



  res <- compute_consensus(b2, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = c(5L, 3L))
  expect_equal(res, structure(list(assessor = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                                3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5), probability = c(0.0333333333333333,
                                                                                                                          0.0333333333333333, 0.0333333333333333, 0.0333333333333333, 0.0333333333333333,
                                                                                                                          0.0333333333333333, 0.0333333333333333, 0.0333333333333333, 0.0333333333333333,
                                                                                                                          0.0333333333333333, 0.0333333333333333, 0.0333333333333333, 0.0333333333333333,
                                                                                                                          0.0333333333333333, 0.0333333333333333, 0.0266666666666667, 0.0266666666666667,
                                                                                                                          0.0266666666666667, 0.0266666666666667, 0.0266666666666667, 0.0266666666666667,
                                                                                                                          0.0266666666666667, 0.0266666666666667, 0.0266666666666667, 0.0266666666666667,
                                                                                                                          0.0266666666666667, 0.0266666666666667, 0.0266666666666667, 0.0266666666666667,
                                                                                                                          0.0266666666666667), item = c("Item 9", "Item 6", "Item 15",
                                                                                                                                                        "Item 10", "Item 3", "Item 1", "Item 7", "Item 11", "Item 13",
                                                                                                                                                        "Item 5", "Item 14", "Item 2", "Item 4", "Item 12", "Item 8",
                                                                                                                                                        "Item 9", "Item 6", "Item 10", "Item 11", "Item 5", "Item 4",
                                                                                                                                                        "Item 3", "Item 1", "Item 15", "Item 8", "Item 2", "Item 7",
                                                                                                                                                        "Item 14", "Item 13", "Item 12"), map_ranking = c(1, 2, 3, 4,
                                                                                                                                                                                                          5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8,
                                                                                                                                                                                                          9, 10, 11, 12, 13, 14, 15)), row.names = c(NA, -30L), class = c("tbl_df",
                                                                                                                                                                                                                                                                          "tbl", "data.frame")))


  res <- compute_consensus(b3, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = c(5L, 3L))

  expect_equal(res, structure(list(assessor = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5,
                                                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5), probability = c(0.03,
                                                                                                        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
                                                                                                        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
                                                                                                        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
                                                                                                        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
                                                ), item = c("Item 9", "Item 9", "Item 6", "Item 6", "Item 15",
                                                            "Item 15", "Item 10", "Item 11", "Item 3", "Item 13", "Item 1",
                                                            "Item 1", "Item 7", "Item 14", "Item 11", "Item 3", "Item 13",
                                                            "Item 7", "Item 5", "Item 10", "Item 14", "Item 5", "Item 2",
                                                            "Item 8", "Item 4", "Item 2", "Item 12", "Item 4", "Item 8",
                                                            "Item 12", "Item 9", "Item 6", "Item 10", "Item 11", "Item 5",
                                                            "Item 4", "Item 3", "Item 1", "Item 15", "Item 8", "Item 2",
                                                            "Item 7", "Item 14", "Item 13", "Item 12"), map_ranking = c(1,
                                                                                                                        1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11,
                                                                                                                        11, 12, 12, 13, 13, 14, 14, 15, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                                                                        10, 11, 12, 13, 14, 15)), row.names = c(NA, -45L), class = c("tbl_df",
                                                                                                                                                                                     "tbl", "data.frame")))
    })
