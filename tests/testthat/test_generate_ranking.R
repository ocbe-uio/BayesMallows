context("Testing generation of initial ranking")

library(dplyr)

# Create some test data
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4
)

pair_comp_tc <- generate_transitive_closure(pair_comp)

test_that("generate_initial_ranking works",{

  expect_error(generate_initial_ranking(pair_comp))
  expect_is(generate_initial_ranking(pair_comp_tc), "matrix")
  expect_true(all(apply(generate_initial_ranking(pair_comp_tc), 1, BayesMallows:::validate_permutation)))
}
)

test_that("generate_initial_ranking with random works",{
  beach_tc <- generate_transitive_closure(beach_preferences)

  small_tc <- beach_tc[beach_tc$assessor %in% 1:6 & beach_tc$bottom_item %in% 1:4 & beach_tc$top_item %in% 1:4, ]
  set.seed(123)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
  expect_equal(init_small, structure(c(3, 1, 1, 1, 4, 3, 1, 4, 3, 4, 2, 1, 4, 3, 2, 3, 3,
                                       4, 2, 2, 4, 2, 1, 2), .Dim = c(6L, 4L)))

  set.seed(321)
  init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
  expect_equal(init_small, structure(c(3, 1, 1, 3, 4, 3, 1, 4, 3, 4, 1, 1, 2, 2, 2, 2, 3,
                                       4, 4, 3, 4, 1, 2, 2), .Dim = c(6L, 4L)))

  expect_error(generate_initial_ranking(pair_comp, random = TRUE))
}
)

