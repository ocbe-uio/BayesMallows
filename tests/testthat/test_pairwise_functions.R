context("Testing functions for pairwise preferences")

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

pair_comp_tc <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1L, 2L,
  1, 2L, 5L,
  1, 4L, 5L,
  1, 1L, 5L,
  2, 1L, 2L,
  2, 2L, 3L,
  2, 3L, 4L,
  2, 1L, 3L,
  2, 1L, 4L,
  2, 2L, 4L
) %>%
  arrange(assessor, bottom_item, top_item)


test_that("transitive closure generation works",{

  pair_comp_returned <- generate_transitive_closure(pair_comp) %>%
    arrange(assessor, bottom_item, top_item)

  expect_equal(pair_comp_tc, pair_comp_returned)

}
)


test_that("Generation of initial ranks works",{
  mat <- generate_initial_ranking(pair_comp_tc)
  expect_equal(TRUE,
               BayesMallows:::validate_initial_ranking(pair_comp_tc, mat)
  )
}
)
