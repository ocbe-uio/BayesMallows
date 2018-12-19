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

