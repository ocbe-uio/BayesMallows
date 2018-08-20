context("Testing misc C++ functions")

test_that(
  "C++ factorial function is correct",
  expect_equal(
    lapply(1:10, BayesMallows:::factorial),
    lapply(1:10, factorial)
  )
)

test_that(
  "C++ binomial coefficient is correct",
  {
    vals <- list(
      c(4, 2),
      c(5, 1),
      c(9, 4),
      c(10, 3)
    )
    expect_equal(
      lapply(X = vals, function(X) BayesMallows:::binomial_coefficient(X[[1]], X[[2]])),
      lapply(X = vals, function(X) choose(X[[1]], X[[2]]))
    )
  }
)

test_that(
  "C++ setdiff is correct",
  {
    x <- 1:100
    y <- sample(x, size = 20)
    expect_equal(
      c(BayesMallows:::std_setdiff(x, y)),
      setdiff(x, y)
    )
  }
)
