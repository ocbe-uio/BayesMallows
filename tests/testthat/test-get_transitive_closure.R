test_that("get_transitive_closure works", {
  dat <- setup_rank_data(preferences = beach_preferences[1:100, ])
  expect_equal(
    as.numeric(get_transitive_closure(dat)[10, 2:3]),
    c(14, 3))

  dd <- data.frame(assessor = 1, bottom_item = 1:2, top_item = 2:1)
  dat <- setup_rank_data(preferences = dd)
  expect_message(
    get_transitive_closure(dat),
    "Intransitive comparisons, no closure exists.")
})
