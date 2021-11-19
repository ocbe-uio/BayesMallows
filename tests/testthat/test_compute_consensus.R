context("Testing compute_consensus")
library(dplyr)
beach_small <- beach_preferences %>%
  filter(bottom_item %in% 1:3, top_item %in% 1:3)

b <- compute_mallows(preferences = beach_preferences, nmc = 100, seed = 123L)
b$burnin <- 2
cp <- compute_consensus(b)
map <- compute_consensus(b, type = "MAP")

test_that("compute_consensus returns correct object", {
  expect_true(inherits(cp, "data.frame"))
  expect_true(inherits(map, "data.frame"))

})

test_that("compute_consensus computes correctly for rho", {
  expect_equal(cp, structure(list(ranking = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                              12, 13, 14, 15), item = c("Item 3", "Item 11", "Item 1", "Item 12",
                                                                        "Item 7", "Item 15", "Item 14", "Item 10", "Item 4", "Item 13",
                                                                        "Item 6", "Item 9", "Item 5", "Item 8", "Item 2"), cumprob = c(0.673469387755102,
                                                                                                                                       0.571428571428571, 0.489795918367347, 0.448979591836735, 0.755102040816327,
                                                                                                                                       0.214285714285714, 0.5, 0.785714285714286, 0.683673469387755,
                                                                                                                                       0.540816326530612, 0.683673469387755, 0.683673469387755, 0.73469387755102,
                                                                                                                                       0.704081632653061, 1)), row.names = c(NA, -15L), class = c("tbl_df",
                                                                                                                                                                                                  "tbl", "data.frame"))
  )

  expect_equal(map,
               structure(list(probability = c(0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653, 0.0510204081632653, 0.0510204081632653,
                                              0.0510204081632653, 0.0510204081632653), item = c("Item 1", "Item 3",
                                                                                                "Item 3", "Item 3", "Item 3", "Item 11", "Item 11", "Item 11",
                                                                                                "Item 1", "Item 7", "Item 11", "Item 12", "Item 1", "Item 7",
                                                                                                "Item 10", "Item 12", "Item 7", "Item 7", "Item 12", "Item 15",
                                                                                                "Item 1", "Item 12", "Item 14", "Item 15", "Item 14", "Item 14",
                                                                                                "Item 15", "Item 15", "Item 4", "Item 9", "Item 10", "Item 10",
                                                                                                "Item 4", "Item 4", "Item 10", "Item 13", "Item 6", "Item 6",
                                                                                                "Item 13", "Item 14", "Item 4", "Item 8", "Item 13", "Item 13",
                                                                                                "Item 6", "Item 6", "Item 9", "Item 9", "Item 5", "Item 5", "Item 8",
                                                                                                "Item 8", "Item 2", "Item 5", "Item 5", "Item 9", "Item 2", "Item 2",
                                                                                                "Item 2", "Item 8"), map_ranking = c(1, 1, 1, 1, 2, 2, 2, 2,
                                                                                                                                     3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8,
                                                                                                                                     8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12,
                                                                                                                                     12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15)), row.names = c(NA,
                                                                                                                                                                                                             -60L), reshapeLong = list(varying = structure(list(map_ranking = c("Item 1",
                                                                                                                                                                                                                                                                                "Item 2", "Item 3", "Item 4", "Item 5", "Item 6", "Item 7", "Item 8",
                                                                                                                                                                                                                                                                                "Item 9", "Item 10", "Item 11", "Item 12", "Item 13", "Item 14",
                                                                                                                                                                                                                                                                                "Item 15")), v.names = "map_ranking", times = c("Item 1", "Item 2",
                                                                                                                                                                                                                                                                                                                                "Item 3", "Item 4", "Item 5", "Item 6", "Item 7", "Item 8", "Item 9",
                                                                                                                                                                                                                                                                                                                                "Item 10", "Item 11", "Item 12", "Item 13", "Item 14", "Item 15"
                                                                                                                                                                                                                                                                                )), v.names = "map_ranking", idvar = c("cluster", "probability",
                                                                                                                                                                                                                                                                                                                       "id"), timevar = "item"), class = "data.frame")
               )
})


b2 <- compute_mallows(preferences = beach_small, nmc = 500, save_aug = TRUE, seed = 123L)
b3 <- compute_mallows(preferences = beach_small, nmc = 500, save_aug = TRUE, aug_thinning = 3, seed = 123L)

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
  expect_equal(res, structure(list(ranking = c(1, 2, 3), item = c("Item 1", "Item 3",
                                                                  "Item 2"), cumprob = c(0.62, 1, 1)), row.names = c(NA, -3L), class = c("tbl_df",
                                                                                                                                         "tbl", "data.frame"))
  )

  res <- compute_consensus(b3, type = "CP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(ranking = c(1, 2, 3), item = c("Item 1", "Item 3",
                                                                  "Item 2"), cumprob = c(0.61, 1, 1)), row.names = c(NA, -3L), class = c("tbl_df",
                                                                                                                                         "tbl", "data.frame"))
  )
  res <- compute_consensus(b2, type = "CP", burnin = 200, parameter = "Rtilde", assessors = c(3L, 5L))

  expect_equal(res, structure(list(assessor = c("3", "3", "3", "5", "5", "5"), ranking = c(1,
                                                                                           2, 3, 1, 2, 3), item = c("Item 1", "Item 3", "Item 2", "Item 1",
                                                                                                                    "Item 3", "Item 2"), cumprob = c(0.696666666666667, 0.836666666666667,
                                                                                                                                                     1, 1, 1, 1)), row.names = c(NA, -6L), class = c("tbl_df", "tbl",
                                                                                                                                                                                                     "data.frame"))
  )

  res <- compute_consensus(b3, type = "CP", burnin = 200, parameter = "Rtilde", assessors = c(3L, 5L))
  expect_equal(res,
               structure(list(assessor = c("3", "3", "3", "5", "5", "5"), ranking = c(1,
                                                                                      2, 3, 1, 2, 3), item = c("Item 1", "Item 3", "Item 2", "Item 1",
                                                                                                               "Item 3", "Item 2"), cumprob = c(0.69, 0.84, 1, 1, 1, 1)), row.names = c(NA,
                                                                                                                                                                                        -6L), class = c("tbl_df", "tbl", "data.frame"))
  )

  res <- compute_consensus(b2, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res, structure(list(probability = c(0.62, 0.62, 0.62), item = c("Item 1",
                                                                               "Item 3", "Item 2"), map_ranking = c(1, 2, 3)), row.names = c(NA,
                                                                                                                                             -3L), reshapeLong = list(varying = structure(list(map_ranking = c("Item 1",
                                                                                                                                                                                                               "Item 2", "Item 3")), v.names = "map_ranking", times = c("Item 1",
                                                                                                                                                                                                                                                                        "Item 2", "Item 3")), v.names = "map_ranking", idvar = c("cluster",
                                                                                                                                                                                                                                                                                                                                 "probability", "id"), timevar = "item"), class = "data.frame"))

  res <- compute_consensus(b3, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = 1L)
  expect_equal(res,
               structure(list(probability = c(0.61, 0.61, 0.61), item = c("Item 1",
                                                                          "Item 3", "Item 2"), map_ranking = c(1, 2, 3)), row.names = c(NA,
                                                                                                                                        -3L), reshapeLong = list(varying = structure(list(map_ranking = c("Item 1",
                                                                                                                                                                                                          "Item 2", "Item 3")), v.names = "map_ranking", times = c("Item 1",
                                                                                                                                                                                                                                                                   "Item 2", "Item 3")), v.names = "map_ranking", idvar = c("cluster",
                                                                                                                                                                                                                                                                                                                            "probability", "id"), timevar = "item"), class = "data.frame")
               )


  res <- compute_consensus(b2, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = c(5L, 3L))
  expect_equal(res, structure(list(assessor = c(3, 3, 3, 5, 5, 5), probability = c(0.533333333333333,
                                                                                   0.533333333333333, 0.533333333333333, 1, 1, 1), item = c("Item 1",
                                                                                                                                            "Item 3", "Item 2", "Item 1", "Item 3", "Item 2"), map_ranking = c(1,
                                                                                                                                                                                                               2, 3, 1, 2, 3)), row.names = c(NA, -6L), reshapeLong = list(varying = structure(list(
                                                                                                                                                                                                                 map_ranking = c("Item 1", "Item 2", "Item 3")), v.names = "map_ranking", times = c("Item 1",
                                                                                                                                                                                                                                                                                                    "Item 2", "Item 3")), v.names = "map_ranking", idvar = c("cluster",
                                                                                                                                                                                                                                                                                                                                                             "probability", "id"), timevar = "item"), class = "data.frame")
               )

  res <- compute_consensus(b3, type = "MAP", burnin = 200, parameter = "Rtilde", assessors = c(5L, 3L))

  expect_equal(
    res,
    structure(list(assessor = c(3, 3, 3, 5, 5, 5), probability = c(0.53,
                                                                   0.53, 0.53, 1, 1, 1), item = c("Item 1", "Item 3", "Item 2",
                                                                                                  "Item 1", "Item 3", "Item 2"), map_ranking = c(1, 2, 3, 1, 2,
                                                                                                                                                 3)), row.names = c(NA, -6L), reshapeLong = list(varying = structure(list(
                                                                                                                                                   map_ranking = c("Item 1", "Item 2", "Item 3")), v.names = "map_ranking", times = c("Item 1",
                                                                                                                                                                                                                                      "Item 2", "Item 3")), v.names = "map_ranking", idvar = c("cluster",
                                                                                                                                                                                                                                                                                               "probability", "id"), timevar = "item"), class = "data.frame")

               )



    })
