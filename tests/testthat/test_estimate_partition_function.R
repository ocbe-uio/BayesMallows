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
#
# test_that(
#   "estimate_partition_function importance sampling gives correct values",
#   {
#     dists <- c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")
#
#     res <- lapply(dists, function(x) {
#                     estimate_partition_function(
#                       method = "importance_sampling",
#                       alpha_vector = seq(from = 1, to = 2, by = .1),
#                       n_items = 10, metric = x,
#                       nmc = 134, degree = 5, seed = 1234)
#                   })
#
#     names(res) <- dists
#
#     expect_equal(
#       res[["footrule"]],
#       c(`(Intercept)` = -30.4691447159505, `I(alpha^1)` = 152.801763425534,
#         `I(alpha^2)` = -210.520483536755, `I(alpha^3)` = 140.114987210433,
#         `I(alpha^4)` = -45.8162389489166, `I(alpha^5)` = 5.89354153793292
#       )
#     )
#
#     expect_equal(
#       res[["spearman"]],
#       c(`(Intercept)` = -73.191387390599, `I(alpha^1)` = 311.594834618677,
#         `I(alpha^2)` = -461.909065750277, `I(alpha^3)` = 331.433424400889,
#         `I(alpha^4)` = -116.214224357447, `I(alpha^5)` = 15.9689321492025
#       )
#     )
#
#     expect_equal(
#       res[["cayley"]],
#       c(`(Intercept)` = 1427.12681864771, `I(alpha^1)` = -4831.94326565112,
#         `I(alpha^2)` = 6513.72488760183, `I(alpha^3)` = -4332.09683073516,
#         `I(alpha^4)` = 1421.9018031319, `I(alpha^5)` = -184.40059946578
#       )
#     )
#
#     expect_equal(
#       res[["hamming"]],
#       c(`(Intercept)` = 1424.14818095233, `I(alpha^1)` = -4820.53503815227,
#         `I(alpha^2)` = 6498.54784270852, `I(alpha^3)` = -4323.86310864666,
#         `I(alpha^4)` = 1420.28219825619, `I(alpha^5)` = -184.386448094199
#       )
#     )
#
#     expect_equal(
#       res[["kendall"]],
#       c(`(Intercept)` = 1162.29711110988, `I(alpha^1)` = -3924.06688047639,
#         `I(alpha^2)` = 5284.33713433091, `I(alpha^3)` = -3512.66960643192,
#         `I(alpha^4)` = 1152.89347081507, `I(alpha^5)` = -149.586139299855
#       )
#     )
#
#     expect_equal(
#       res[["ulam"]],
#       c(`(Intercept)` = 1353.64045264671, `I(alpha^1)` = -4575.69490261372,
#         `I(alpha^2)` = 6163.32301333952, `I(alpha^3)` = -4096.48687926769,
#         `I(alpha^4)` = 1343.95069222954, `I(alpha^5)` = -174.238719114932
#       )
#     )
#
#       }
# )
