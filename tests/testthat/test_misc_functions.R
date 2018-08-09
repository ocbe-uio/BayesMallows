context("Testing misc functions")

test_that(
  "validate_permutation is correct",{
    expect_equal(
      BayesMallows:::validate_permutation(c(1, 3, 3)),
      FALSE
    )
    expect_equal(
      BayesMallows:::validate_permutation(c(1, 2, 3)),
      TRUE
    )
    expect_equal(
      BayesMallows:::validate_permutation(c(1, NA_real_, 3, NA_real_)),
      TRUE
    )
    expect_equal(
      BayesMallows:::validate_permutation(c(NA_real_, NA_real_, NA_real_)),
      TRUE
    )
  }

)
