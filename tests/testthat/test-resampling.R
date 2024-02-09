test_that("resampling options work", {
  set.seed(1)
  mod0 <- sample_prior(1000, 20)

  for (r in c("stratified", "systematic", "residual", "multinomial")) {
    mod1 <- update_mallows(
      model = mod0,
      new_data = setup_rank_data(potato_weighing[1:3, ]),
      smc_options = set_smc_options(n_particles = 20, resampler = r)
    )
    expect_s3_class(mod1, "SMCMallows")
  }

  expect_error(
    update_mallows(
      model = mod0,
      new_data = setup_rank_data(potato_weighing[1:3, ]),
      smc_options = set_smc_options(n_particles = 20, resampler = "gaussian")
    ),
    "'arg' should be one of"
  )
})
