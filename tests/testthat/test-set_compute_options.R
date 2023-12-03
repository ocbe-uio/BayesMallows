test_that("set_compute_options works", {
  s <- set_compute_options()
  expect_s3_class(s, "BayesMallowsComputeOptions")
  expect_error(set_compute_options(nmc = -1))
  expect_error(set_compute_options(burnin = -1))
  expect_error(set_compute_options(nmc = 3, burnin = 5))
  expect_error(set_compute_options(alpha_prop_sd = -1))
  expect_error(set_compute_options(aug_thinning = 1e6))
  expect_error(set_compute_options(clus_thinning = 1e6))
  expect_error(set_compute_options(rho_thinning = 1e6))
  expect_error(set_compute_options(alpha_jump = 1e6))

  f <- file()
  write("yes", f)
  options("ask_opts.con" = f)
  s <- set_compute_options(save_ind_clus = TRUE)
  close(f)
  expect_equal(s$save_ind_clus, TRUE)

  f <- file()
  write("no", f)
  options("ask_opts.con" = f)
  expect_error(
    set_compute_options(save_ind_clus = TRUE),
    "quitting")
  close(f)
})
