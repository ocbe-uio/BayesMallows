test_that("set_model_options works", {
  expect_error(
    set_model_options(metric = "manhattan"),
    "'arg' should be one of")

  expect_error(
    set_model_options(error_model = "sturm liouville"),
    "'arg' should be one of")

  expect_s3_class(set_model_options(), "BayesMallowsModelOptions")
})
