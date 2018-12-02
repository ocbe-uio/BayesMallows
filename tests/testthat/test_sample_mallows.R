context("Testing sample_mallows")

set.seed(1)
# Number of items
n_items <- 15
# Set the consensus ranking
rho0 <- seq(from = 1, to = n_items, by = 1)
# Set the scale
alpha0 <- 10

for(m in c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")){
  samples <- sample_mallows(rho0 = rho0, alpha0 = alpha0, n_samples = 100,
                            burnin = 1000, thinning = 1000, metric = m)
  test_that(
    "sample_mallows returns correct values",
    expect_true(mean(samples[, 1]) < mean(samples[, n_items]))
  )

  test_that(
    "sample_mallows returns matrix",
    expect_is(samples, "matrix")
  )
}

# Check that an improper ranking yields error message
test_that(
  "sample_mallows handles wrong input correctly", {
    expect_error(sample_mallows(rho0 = c(1, 1, 2), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2, NA_real_), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = 1, diagnostic = TRUE))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = -2))
  }
)

