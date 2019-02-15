library(dplyr)
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


test_that("compute_mallows error model works", {
  preferences <- data.frame(assessor = c(1, 1, 2, 2),
                            bottom_item = c(1, 2, 1, 2),
                            top_item = c(2, 1, 2, 3)
                            )
  expect_error(invisible(capture.output(compute_mallows(preferences = preferences, nmc = 10))))
  expect_s3_class(compute_mallows(preferences = preferences, error_model = "bernoulli", nmc = 10),
                  "BayesMallows")

})

test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(rankings = mat, nmc = 30)
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")

})


test_that("compute_mallows runs with the right distances", {
  for(metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")){
    expect_s3_class(compute_mallows(potato_visual, metric = metric, nmc = 3), "BayesMallows")
  }

})

test_that("compute_mallows handles integer preferences", {
  m <- beach_preferences %>%
    mutate_all(as.integer) %>%
    compute_mallows(preferences = .)

})

test_that("compute_mallows handles data with lots of missings",{
  R_partial2 <- structure(c(NA, NA, NA, NA, NA, NA, 9, NA, NA, 7, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, 7, 8, 10, NA, NA, 9, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6, NA, 5, 2,
                            6, 5, 6, 6, 5, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, 3, 4, NA, 3, 3, 3, 3, 4, 5, 3, 3, 3, NA, 3, 3, 4, NA,
                            7, 8, 3, 3, 10, 5, 4, NA, NA, NA, 8, NA, NA, NA, NA, NA, 11,
                            NA, NA, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 4, 2,
                            2, 2, 4, 2, 2, 2, NA, NA, 4, 7, 5, 4, 6, 7, 2, 6, 6, 7, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9,
                            7, 8, NA, 8, 10, 6, NA, 5, NA, 6, 6, 5, 4, 5, NA, 4, 4, 5, NA,
                            NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, 9, NA), .Dim = c(12L, 20L))

  m <- compute_mallows(R_partial2)
  expect_s3_class(assess_convergence(m), "gg")

}
          )
