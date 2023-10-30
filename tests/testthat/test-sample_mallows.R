set.seed(1)
# Number of items
n_items <- 15
# Set the consensus ranking
rho0 <- seq(from = 1, to = n_items, by = 1)
# Set the scale
alpha0 <- 10

for (m in c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")) {
  samples <- sample_mallows(
    rho0 = rho0, alpha0 = alpha0, n_samples = 100,
    burnin = 1000, thinning = 1000, metric = m, leap_size = 1
  )
  test_that(
    "sample_mallows returns correct values",
    {
      expect_true(mean(samples[, 1]) < mean(samples[, n_items]))
    }

  )

  test_that(
    "sample_mallows returns matrix",
    {
      expect_true(inherits(samples, "matrix"))
    }
  )
}

# Check that an improper ranking yields error message
test_that(
  "sample_mallows handles wrong input correctly",
  {
    expect_error(sample_mallows(rho0 = c(1, 1, 2), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2, NA_real_), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = 1, diagnostic = TRUE))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = -2))
  }
)

test_that("sample_mallows diagnostics work", {
  set.seed(1)
  expect_message(
    {
      pdf(NULL)
      test <- sample_mallows(
        rho0 = rho0, alpha0 = alpha0, diagnostic = TRUE,
        n_samples = 1000, burnin = 1, thinning = 1
      )
    },
    "Items not provided by user. Picking 5 at random."
  )
  expect_equal(
    test[c(200, 300, 400), ],
    structure(c(
      1, 1, 1, 2, 3, 2, 3, 2, 4, 6, 4, 3, 5, 5, 5, 4, 6,
      6, 7, 8, 9, 11, 9, 10, 9, 7, 7, 8, 10, 12, 12, 11, 11, 10, 12,
      8, 13, 14, 13, 14, 13, 14, 15, 15, 15
    ), .Dim = c(3L, 15L), .Dimnames = list(
      NULL, c(
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15"
      )
    ))
  )
})
