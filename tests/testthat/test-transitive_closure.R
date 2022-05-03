context("Testing functions for pairwise preferences")

# Create some test data
pair_comp <- data.frame(
  assessor = c(1, 1, 1, 2, 2, 2),
  bottom_item = c(1, 2, 4, 1, 2, 3),
  top_item = c(2, 5, 5, 2, 3, 4)
)

pair_comp_tc <- data.frame(
  assessor = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
  bottom_item = c(1, 1, 2, 4, 1, 1, 1, 2, 2, 3),
  top_item = c(2, 5, 5, 5, 2, 3, 4, 3, 4, 4)
)

class(pair_comp_tc) <- c("BayesMallowsTC", class(pair_comp_tc))

test_that("transitive closure generation works", {
  pair_comp_returned <- generate_transitive_closure(pair_comp)
  pair_comp_returned <- pair_comp_returned[order(
    pair_comp_returned$assessor,
    pair_comp_returned$bottom_item,
    pair_comp_returned$top_item
  ), ]

  expect_equal(
    as.matrix(pair_comp_returned, rownames.force = FALSE),
    as.matrix(pair_comp_tc, rownames.force = FALSE)
  )
})

test_that("transitive closure generation discovers inconsistencies", {
  pair_comp_inc <- rbind(
    pair_comp,
    data.frame(assessor = 1, bottom_item = 5L, top_item = 2L)
  )
  expect_error(invisible(capture.output(generate_transitive_closure(pair_comp_inc))))
})
