set.seed(1)
# Number of items
n_items <- 15
# Set the consensus ranking
rho0 <- seq(from = 1, to = n_items, by = 1)
# Set the scale
alpha0 <- 10

test_that("sample_mallows works with all distances", {
  for (m in c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")) {
    print(m)
    samples <- sample_mallows(
      rho0 = rho0, alpha0 = alpha0, n_samples = 200,
      burnin = 1000, thinning = 100, metric = m, leap_size = 1
    )
    expect_true(mean(samples[, 1]) < mean(samples[, n_items]))
    expect_true(inherits(samples, "matrix"))

    f <- file()
    write("\n", f)
    options("ask_opts.con" = f)
    expect_output(
      sample_mallows(
        rho0 = rho0, alpha0 = alpha0, n_samples = 1000,
        metric = m, diagnostic = TRUE,
        items_to_plot = 1:4
      ),
      "to see the next plot")
    close(f)
  }
})

test_that(
  "sample_mallows handles wrong input correctly",
  {
    expect_error(sample_mallows(rho0 = c(1, 1, 2), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2, NA_real_), alpha0 = 1, n_samples = 1))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = 1, diagnostic = TRUE))
    expect_error(sample_mallows(rho0 = c(1, 2), alpha0 = 1, n_samples = -2))
  }
)
