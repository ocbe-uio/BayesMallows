test_that("expected dist works", {
  expect_equal(round(expected_dist(1,5,metric="kendall"), 6), 1.749137)
  expect_equal(round(expected_dist(2,6,metric="cayley"), 6), 1.375779)
  expect_equal(round(expected_dist(1.5,7,metric="hamming"), 6), 2.69246)
  expect_equal(round(expected_dist(5,30,"ulam"), 6), 4.133538)
  expect_equal(round(expected_dist(3.5,45,"footrule"), 6), 0.080459)
  expect_equal(round(expected_dist(4,10,"spearman"), 6), 0.006033)
})
