test_that("progress reporter works", {
  expect_output(
    mm <- compute_mallows(
      data = setup_rank_data(potato_visual),
      compute_options = set_compute_options(nmc = 4),
      progress_report = set_progress_report(verbose = TRUE, report_interval = 2)
    ),
    "First 2 iterations of Metropolis-Hastings algorithm completed."
  )

  expect_error(
    set_progress_report(1),
    "verbose must be a logical value of length one"
  )

  expect_error(
    set_progress_report(TRUE, -3),
    "report_interval must be a strictly positive number of length one"
  )
})
