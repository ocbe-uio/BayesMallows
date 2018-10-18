context("Testing compute_mallows")

test_that("miscellaneous input validation", {
  expect_error(compute_mallows(nmc = 1000, alpha_prop_sd = 1))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, alpha_jump = 102))
  expect_error(compute_mallows(rankings = potato_visual, lambda = 0))
  expect_error(compute_mallows(rankings = potato_visual, lambda = -10))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, rho_thinning = 200))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, aug_thinning = 200))
  expect_error(compute_mallows(rankings = potato_visual, nmc = -100))
})

test_that("rho_init is properly validated",{
  m <- potato_visual
  expect_error(compute_mallows(rankings = m, rho_init = 1:(ncol(m) - 1)))
  expect_error(compute_mallows(rankings = m, rho_init = c(potato_true_ranking[-1], 22)))
  expect_error(compute_mallows(rankings = m, rho_init = c(NA_real_, 2:ncol(m))))
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
