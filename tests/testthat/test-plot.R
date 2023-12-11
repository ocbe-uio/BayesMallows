test_that("plot.SMCMallows works", {
  set.seed(1)
  data_first_batch <- potato_visual[1:4, ]

  mod_init <- compute_mallows(
    data = setup_rank_data(data_first_batch),
    compute_options = set_compute_options(nmc = 100, burnin = 0)
  )

  data_second_batch <- potato_visual[5:8, ]
  mod_second <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = data_second_batch),
    smc_options = set_smc_options(n_particles = 10)
  )

  p <- plot(mod_second)
  expect_s3_class(p, "ggplot")
  expect_equal(dim(p$data), c(10, 4))
  p <- plot(mod_second, parameter = "rho", items = 1:4)
  expect_equal(dim(p$data), c(16, 5))
  expect_message(
    p <- plot(mod_second, parameter = "rho"),
    "Items not provided by user. Picking 5 at random."
  )

  data_third_batch <- potato_visual[9:12, ]
  mod_final <- update_mallows(
    model = mod_second, new_data = setup_rank_data(rankings = data_third_batch)
  )

  p <- plot(mod_final)
  expect_equal(p$labels$y, "Posterior density")
  expect_equal(dim(p$data), c(10, 4))

  p <- plot(mod_final, parameter = "rho", items = c("P19", "P8"))
  expect_s3_class(p, "ggplot")
  expect_equal(dim(p$data), c(8, 5))
  expect_equal(as.character(unique(p$data$item)), c("P19", "P8"))

  expect_error(
    plot(mod_final, parameter = "theta"),
    "'arg' should be one of"
  )
})
