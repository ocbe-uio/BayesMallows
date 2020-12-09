context("Testing function estimate_partition_function")


test_that(
  "estimate_partition_function fails for wrong asymptotic metrics",
  {
    for(metric in c("cayley", "hamming", "kendall", "ulam")){
      expect_error(
        estimate_partition_function(method = "asymptotic",
                                    alpha_vector = seq(from = 1, to = 2, by = .1),
                                    n_items = 10,
                                    metric = metric,
                                    n_iterations = 50,
                                    K = 20, degree = 5)
      )
    }

  }
)

test_that(
  "estimate_partition_function asymptotic gives correct values",
  {
    expect_equal(
      estimate_partition_function(method = "asymptotic",
                                  alpha_vector = seq(from = 1, to = 2, by = .1),
                                  n_items = 10, metric = "footrule",
                                  n_iterations = 50, K = 20, degree = 5),
      c(`(Intercept)` = 15.1041186472153, `I(alpha^1)` = -3.32366499344578,
        `I(alpha^2)` = 0.221067735142993, `I(alpha^3)` = 0.00983874213603723,
        `I(alpha^4)` = -0.00400646965806447, `I(alpha^5)` = 0.000290662453972935
      )
    )
    expect_equal(
      estimate_partition_function(method = "asymptotic",
                                  alpha_vector = seq(from = 1, to = 2, by = .2),
                                  n_items = 20, metric = "spearman",
                                  n_iterations = 55, K = 21, degree = 4),
      c(`(Intercept)` = 34.3498359623487, `I(alpha^1)` = -22.0673491727874,
        `I(alpha^2)` = 10.8140117937552, `I(alpha^3)` = -3.18806119268205,
        `I(alpha^4)` = 0.394939624914785)
    )
  }
)

test_that(
  "estimate_partition_function importance sampling gives correct values",
  {
    dists <- c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")
    set.seed(1234)
    res <- lapply(dists, function(x) {
                    estimate_partition_function(
                      method = "importance_sampling",
                      alpha_vector = seq(from = 1, to = 2, by = .1),
                      n_items = 10, metric = x,
                      nmc = 134, degree = 5)
                  })

    names(res) <- dists

    expect_equal(
      res[["footrule"]],
      c(`(Intercept)` = -30.4691447159505, `I(alpha^1)` = 152.801763425534,
        `I(alpha^2)` = -210.520483536755, `I(alpha^3)` = 140.114987210433,
        `I(alpha^4)` = -45.8162389489166, `I(alpha^5)` = 5.89354153793292
      )
    )

    expect_equal(
      res[["spearman"]],
      c(`(Intercept)` = 2.15976899368893, `I(alpha^1)` = 10.6313257568846,
        `I(alpha^2)` = 4.24514443129334, `I(alpha^3)` = -18.8229864455337,
        `I(alpha^4)` = 11.6710538209137, `I(alpha^5)` = -2.2215884037716
      )
    )

    expect_equal(
      res[["cayley"]],
      c(`(Intercept)` = -326.499454772294, `I(alpha^1)` = 1177.20606429431,
        `I(alpha^2)` = -1630.78697011216, `I(alpha^3)` = 1129.41670099163,
        `I(alpha^4)` = -390.975111187882, `I(alpha^5)` = 53.990569504206
      )
    )

    expect_equal(
      res[["hamming"]],
      c(`(Intercept)` = -333.646767315888, `I(alpha^1)` = 1171.74149116192,
        `I(alpha^2)` = -1555.27700698541, `I(alpha^3)` = 1015.11114448604,
        `I(alpha^4)` = -326.941349359867, `I(alpha^5)` = 41.6154937354296
      )
    )

    expect_equal(
      res[["kendall"]],
      c(`(Intercept)` = -327.821301817314, `I(alpha^1)` = 1098.11871012986,
        `I(alpha^2)` = -1382.93367276787, `I(alpha^3)` = 848.61070203972,
        `I(alpha^4)` = -254.119854587992, `I(alpha^5)` = 29.6918584845644
      )
    )

    expect_equal(
      res[["ulam"]],
      c(`(Intercept)` = 1539.56192689499, `I(alpha^1)` = -5418.3727314399,
        `I(alpha^2)` = 7600.30334960907, `I(alpha^3)` = -5262.97091487374,
        `I(alpha^4)` = 1797.77976372225, `I(alpha^5)` = -242.277562572955
      )
    )

      }
)
