test_that(
  "estimate_partition_function fails for wrong asymptotic metrics",
  {
    for (metric in c("cayley", "hamming", "kendall", "ulam")) {
      expect_error(
        estimate_partition_function(
          method = "asymptotic",
          alpha_vector = seq(from = 1, to = 2, by = .1),
          n_items = 10,
          metric = metric,
          n_iterations = 50,
          K = 20, degree = 5
        )
      )
    }
  }
)

test_that(
  "estimate_partition_function asymptotic gives correct values",
  {
    expect_equal(
      estimate_partition_function(
        method = "asymptotic",
        alpha_vector = seq(from = 1, to = 2, by = .1),
        n_items = 10, metric = "footrule",
        n_iterations = 50, K = 20, degree = 5
      ),
      c(
        `(Intercept)` = 15.1041186472153, `I(alpha^1)` = -3.32366499344578,
        `I(alpha^2)` = 0.221067735142993, `I(alpha^3)` = 0.00983874213603723,
        `I(alpha^4)` = -0.00400646965806447, `I(alpha^5)` = 0.000290662453972935
      )
    )
    expect_equal(
      estimate_partition_function(
        method = "asymptotic",
        alpha_vector = seq(from = 1, to = 2, by = .2),
        n_items = 20, metric = "spearman",
        n_iterations = 55, K = 21, degree = 4
      ),
      c(
        `(Intercept)` = 34.3498359623487, `I(alpha^1)` = -22.0673491727874,
        `I(alpha^2)` = 10.8140117937552, `I(alpha^3)` = -3.18806119268205,
        `I(alpha^4)` = 0.394939624914785
      )
    )
  }
)
