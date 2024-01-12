test_that("set_initial_values works", {
  expect_s3_class(
    set_initial_values(), "BayesMallowsInitialValues"
  )
  expect_error(
    set_initial_values(rho_init = c(1, 1, 2)),
    "rho_init must be a proper permutation"
  )
  expect_error(
    set_initial_values(rho_init = c(NA, 1)),
    "rho_init cannot have missing values"
  )
})
