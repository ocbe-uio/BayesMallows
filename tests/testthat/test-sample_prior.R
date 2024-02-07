test_that("sample_prior works", {
  set.seed(1)
  mm <- sample_prior(3, 10)
  expect_s3_class(mm, "BayesMallowsPriorSamples")
  expect_equal(mm$n_items, 10)
  expect_equal(mm$items, 1:10)
  expect_equal(dim(mm$rho), c(10, 3))
  expect_equal(length(mm$alpha), 3)
  expect_equal(mm$alpha[[2]], 1882.40159634668)
  expect_equal(mm$rho[2, 2], 5)
})

test_that("scale and shape matter", {
  set.seed(1)
  mm1 <- sample_prior(100, 20)
  mm2 <- sample_prior(100, 20, priors = set_priors(gamma = 2, lambda = .1))

  expect_lt(mean(mm2$alpha), mean(mm1$alpha))
})
