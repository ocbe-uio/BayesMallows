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
  "C++ sample_int is correct",
  {

    n <- 10000L
    probs <- c(0.1, 0.2, 0.7)
    values <- integer(n)
    for(i in seq(1L, n, 1L)){
      values[[i]] <- BayesMallows:::sample_int(probs)
    }
    freqs <- table(values) / length(values)

    expect_equal(sum(freqs), 1)

    # Skip the following on CRAN, since they may not always be
    # TRUE, due to the randomness in sampling
    skip_on_cran()
    diff <- abs(probs - freqs)
    for(i in 1:3){
      expect_lt(diff[[i]], 0.02)
    }


  }
)
