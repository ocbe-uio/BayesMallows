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

test_that(
  "validate_rank_ordering is correct", {

    # Here is a rank matrix
    rankings <- matrix(c(5, 5, 4,
                         2, 4, 1,
                         4, 3, 5,
                         3, 2, 2,
                         1, 1, 3),
                       nrow = 3)

    # Here is a linear ordering that is consistent with the matrix
    linear_ordering <- list(
      c(1, 2),
      c(1, 2, 3, 4, 5),
      numeric()
    )

    expect_true(
      validate_rank_ordering(rankings, linear_ordering)
    )

    # Here is another linear ordering that is consistent with the matrix
    linear_ordering <- list(
      numeric(),
      numeric(),
      numeric()
    )

    expect_true(
      validate_rank_ordering(rankings, linear_ordering)
    )

    # Here is another linear ordering that is consistent with the matrix
    linear_ordering <- list(
      c(1, 3, 4, 2, 5),
      c(1, 4, 5),
      c(3, 1, 5, 4, 2)
    )

    expect_true(
      validate_rank_ordering(rankings, linear_ordering)
    )

    # Here is a linear ordering that is inconsistent with the matrix
    linear_ordering <- list(
      c(5, 3, 4, 2, 1), # This one is inconsistent
      c(1, 4, 5),
      c(3, 1, 5, 4, 2)
    )

    expect_false(
      validate_rank_ordering(rankings, linear_ordering)
    )

    # Here is another linear ordering that is inconsistent with the matrix
    linear_ordering <- list(
      c(5, 3, 4, 2, 1), # This one is inconsistent
      c(5, 1), # This one is inconsistent
      c(3, 1, 5, 4, 2)
    )

    expect_false(
      validate_rank_ordering(rankings, linear_ordering)
    )

  }
)
