test_that("expected dist works", {
  expect_equal(round(expected_dist(5,5,metric="kendall"), 6), 1.749137)
  expect_equal(round(expected_dist(12,6,metric="cayley"), 6), 1.375779)
  expect_equal(round(expected_dist(1.5 * 7,7,metric="hamming"), 6), 2.69246)
  expect_equal(round(expected_dist(5 * 30,30,"ulam"), 6), 4.133538)
  expect_equal(round(expected_dist(3.5 * 45,45,"footrule"), 6), 0.080459)
  expect_equal(round(expected_dist(4 * 10,10,"spearman"), 6), 0.006033)
})


test_that("expected dist fails when it should", {
  expect_error(expected_dist(10, 15, "spearman"))
  expect_error(expected_dist(10, 150, "footrule"))
  expect_error(expected_dist(10, 150, "ulam"))
})
