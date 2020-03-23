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
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  1, 1, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4,
  2, 1, 3,
  2, 1, 4,
  2, 2, 4
) %>%
  arrange(assessor, bottom_item, top_item)

class(pair_comp_tc) <- c("BayesMallowsTC", class(pair_comp_tc))


test_that("transitive closure generation works",{

  pair_comp_returned <- generate_transitive_closure(pair_comp) %>%
    arrange(assessor, bottom_item, top_item)

  expect_equal(as.matrix(pair_comp_tc), as.matrix(pair_comp_returned))

}
)

test_that("transitive closure generation discovers inconsistencies",{
  pair_comp_inc <- bind_rows(pair_comp,
                             tibble(assessor = 1, bottom_item = 5L, top_item = 2L))
  expect_error(invisible(capture.output(generate_transitive_closure(pair_comp_inc))))
})


