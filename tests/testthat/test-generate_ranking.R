context("Testing generation of initial ranking")

# Create some test data
pair_comp <- data.frame(
  assessor = c(rep(1, 3), rep(2, 3)),
  bottom_item = c(1, 2, 4, 1, 2, 3),
  top_item = c(2, 5, 5, 2, 3, 4)
)

pair_comp_tc <- generate_transitive_closure(pair_comp)
beach_tc <- generate_transitive_closure(beach_preferences)

test_that("generate_initial_ranking works", {
  expect_error(generate_initial_ranking(pair_comp))
  expect_is(generate_initial_ranking(pair_comp_tc), "matrix")
  expect_true(all(apply(generate_initial_ranking(pair_comp_tc), 1, BayesMallows:::validate_permutation)))
  expect_error(generate_initial_ranking(beach_tc, n_items = 10))
})

test_that("generate_initial_ranking with shuffle_unranked works", {
  small_tc <- beach_tc[beach_tc$assessor %in% 1:6 & beach_tc$bottom_item %in% 1:4 & beach_tc$top_item %in% 1:4, ]
  set.seed(123)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, shuffle_unranked = TRUE)
  expect_equal(init_small, structure(c(
    2L, 1L, 2L, 2L, 3L, 2L, 3L, 3L, 3L, 4L, 1L, 4L, 1L,
    4L, 1L, 3L, 4L, 1L, 4L, 2L, 4L, 1L, 2L, 3L
  ), .Dim = c(6L, 4L), dimnames = list(
    c("1", "2", "3", "4", "5", "6"), NULL)))

  set.seed(321)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
  expect_equal(init_small, structure(c(
    3, 1, 1, 3, 4, 3, 1, 4, 3, 4, 1, 1, 2, 2, 2, 2, 3,
    4, 4, 3, 4, 1, 2, 2
  ), .Dim = c(6L, 4L), dimnames = list(c("1",
                                         "2", "3", "4", "5", "6"), NULL)))

  expect_error(generate_initial_ranking(pair_comp, shuffle_unranked = TRUE))
})



test_that("generate_initial_ranking with random works", {
  small_tc <- beach_tc[beach_tc$assessor %in% 1:6 & beach_tc$bottom_item %in% 1:4 & beach_tc$top_item %in% 1:4, ]
  set.seed(123)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
  expect_equal(init_small, structure(c(
    3, 1, 1, 1, 4, 3, 1, 4, 3, 4, 2, 1, 4, 3, 2, 3, 3,
    4, 2, 2, 4, 2, 1, 2
  ), .Dim = c(6L, 4L), dimnames = list(c("1",
                                         "2", "3", "4", "5", "6"), NULL)))

  set.seed(321)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
  expect_equal(init_small, structure(c(
    3, 1, 1, 3, 4, 3, 1, 4, 3, 4, 1, 1, 2, 2, 2, 2, 3,
    4, 4, 3, 4, 1, 2, 2
  ), .Dim = c(6L, 4L),
  dimnames = list(c("1", "2", "3", "4", "5", "6"), NULL)))

  expect_error(generate_initial_ranking(pair_comp, random = TRUE))

  expect_error(generate_initial_ranking(
    tc = small_tc, n_items = 4, random = TRUE,
    random_limit = 3
  ))
  expect_error(generate_initial_ranking(tc = beach_tc, random = TRUE))
})
