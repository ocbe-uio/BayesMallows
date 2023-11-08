# Create some test data
pair_comp <- data.frame(
  assessor = c(rep(1, 3), rep(2, 3)),
  bottom_item = c(1, 2, 4, 1, 2, 3),
  top_item = c(2, 5, 5, 2, 3, 4)
)

dat <- setup_rank_data(preferences = pair_comp)

test_that("generate_initial_ranking works", {
  expect_type(dat$rankings, "integer")
  expect_true(all(apply(dat$rankings, 1, validate_permutation)))
})

beach_small <- subset(beach_preferences, assessor %in% 1:6 & bottom_item %in% 1:4 & top_item %in% 1:4)
test_that("generate_initial_ranking with shuffle_unranked works", {
  set.seed(123)
  dat <- setup_rank_data(preferences = beach_small, shuffle_unranked = TRUE, random = FALSE)
  expect_equal(dat$rankings, structure(c(
    1L, 4L, 1L, 1L, 2L, 4L, 2L, 4L, 2L, 4L, 3L, 1L, 2L,
    4L, 1L, 2L, 3L, 3L, 3L, 3L
  ), dim = 5:4, dimnames = list(c(
    "2",
    "3", "4", "5", "6"
  ), NULL)))

  set.seed(321)
  dat <- setup_rank_data(preferences = beach_small, shuffle_unranked = TRUE, random = FALSE)
  expect_equal(dat$rankings, structure(c(
    1L, 1L, 3L, 4L, 2L, 4L, 3L, 4L, 1L, 4L, 2L, 2L, 2L,
    3L, 1L, 3L, 4L, 1L, 2L, 3L
  ), dim = 5:4, dimnames = list(c(
    "2",
    "3", "4", "5", "6"
  ), NULL)))
})



test_that("generate_initial_ranking with random works", {
  set.seed(123)
  dat <- setup_rank_data(preferences = beach_small, random = TRUE)
  expect_equal(dat$rankings, structure(c(
    1L, 3L, 1L, 1L, 1L, 4L, 2L, 4L, 4L, 3L, 3L, 1L, 3L,
    3L, 2L, 2L, 4L, 2L, 2L, 4L
  ), dim = 5:4, dimnames = list(c(
    "2",
    "3", "4", "5", "6"
  ), NULL)))
})
