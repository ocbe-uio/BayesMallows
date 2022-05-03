test_that("compute posterior intervals works", {
  set.seed(3344)
  m <- compute_mallows(potato_visual, nmc = 10)
  expect_error(compute_posterior_intervals(m))
  expect_error(compute_posterior_intervals(m, burnin = 100))
  expect_error(compute_posterior_intervals(m, burnin = 7, parameter = "dsdsd"))
  expect_error(compute_posterior_intervals(m, burnin = 7, level = 1.2))
  expect_equal(
    compute_posterior_intervals(m, burnin = 7, level = .05, parameter = "alpha"),
    structure(list(
      parameter = "alpha", mean = 0.73, median = 0.744,
      conf_level = "5 %", hpdi = "[0.699,0.699]", central_interval = "[0.741,0.744]"
    ), row.names = c(
      NA,
      -1L
    ), class = "data.frame")
  )

  expect_equal(
    compute_posterior_intervals(m, burnin = 7, level = .1, parameter = "cluster_probs"),
    structure(list(
      parameter = "cluster_probs", mean = 1, median = 1,
      conf_level = "10 %", hpdi = "[1.000,1.000]", central_interval = "[1.000]"
    ), row.names = c(
      NA,
      -1L
    ), class = "data.frame")
  )

  expect_equal(
    compute_posterior_intervals(m, burnin = 7, level = .01, parameter = "rho"),
    structure(list(
      item = structure(1:20, .Label = c(
        "P1", "P2",
        "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12",
        "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20"
      ), class = "factor"),
      parameter = c(
        "rho", "rho", "rho", "rho", "rho", "rho", "rho",
        "rho", "rho", "rho", "rho", "rho", "rho", "rho", "rho", "rho",
        "rho", "rho", "rho", "rho"
      ), mean = c(
        18, 15, 5, 20, 10,
        11, 8, 7, 3, 14, 2, 17, 12, 19, 12, 1, 9, 16, 4, 6
      ), median = c(
        18,
        15, 5, 20, 10, 11, 8, 7, 3, 14, 2, 17, 12, 19, 13, 1, 9,
        16, 4, 6
      ), conf_level = c(
        "1 %", "1 %", "1 %", "1 %", "1 %",
        "1 %", "1 %", "1 %", "1 %", "1 %", "1 %", "1 %", "1 %", "1 %",
        "1 %", "1 %", "1 %", "1 %", "1 %", "1 %"
      ), hpdi = c(
        "[18]",
        "[15]", "[5]", "[20]", "[10]", "[11]", "[8]", "[7]", "[3]",
        "[14]", "[2]", "[17]", "[12]", "[19]", "[13]", "[1]", "[9]",
        "[16]", "[4]", "[6]"
      ), central_interval = c(
        "[18]", "[15]",
        "[5]", "[20]", "[10]", "[11,11]", "[8]", "[7]", "[3]", "[14]",
        "[2]", "[17]", "[12,12]", "[19]", "[13,13]", "[1]", "[9]",
        "[16]", "[4]", "[6]"
      )
    ), row.names = c(NA, -20L), class = "data.frame")
  )

  set.seed(22)
  m <- compute_mallows(potato_visual, nmc = 10, n_clusters = 2)
  expect_equal(
    compute_posterior_intervals(m, burnin = 8),
    structure(list(cluster = structure(1:2, .Label = c(
      "Cluster 1",
      "Cluster 2"
    ), class = "factor"), parameter = c("alpha", "alpha"), mean = c(0.57, 0.817), median = c(0.57, 0.817), conf_level = c(
      "95 %",
      "95 %"
    ), hpdi = c("[0.520,0.620]", "[0.751,0.884]"), central_interval = c(
      "[0.522,0.617]",
      "[0.754,0.880]"
    )), class = "data.frame", row.names = c(NA, -2L))
  )
})
