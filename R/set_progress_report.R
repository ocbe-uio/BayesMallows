#' @title Set progress report options for MCMC algorithm
#'
#' @description Specify whether progress should be reported, and how often.
#'
#' @param verbose Boolean specifying whether to report progress or not. Defaults
#'   to `FALSE`.
#'
#' @param report_interval Strictly positive number specifying how many
#'   iterations of MCMC should be run between each progress report. Defaults to
#'   `1000`.
#'
#' @return An object of class `"BayesMallowsProgressReport"`, to be provided in
#'   the `progress_report` argument to [compute_mallows()] and
#'   [compute_mallows_mixtures()].
#' @export
#'
#' @references \insertAllCited{}
#'
#' @family preprocessing
#'
set_progress_report <- function(verbose = FALSE, report_interval = 1000) {
  validate_positive(report_interval)
  validate_logical(verbose)
  ret <- as.list(environment())
  class(ret) <- "BayesMallowsProgressReport"
  ret
}
