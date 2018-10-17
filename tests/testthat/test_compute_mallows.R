context("Testing compute_mallows")

test_that("rho_init is properly validated",{
  m <- potato_visual
  expect_error(compute_mallows(rankings = m, rho_init = 1:(ncol(m) - 1)))
  expect_error(compute_mallows(rankings = m, rho_init = c(potato_true_ranking[-1], 22)))
}
)

test_that("compute_mallows discovers inconsistent rankings",{
    expect_error(compute_mallows(
      rankings = matrix(c(1, 2, -3,
                          1, 2, 3), nrow = 2, byrow = TRUE)
    ))
  expect_error(compute_mallows(
    rankings = matrix(c(1, 2, 3,
                        1, 2, 2), nrow = 2, byrow = TRUE)
  ))
  }
)
