#' @title SMC-Mallows New Users
#' @description Function to perform resample-move SMC algorithm where we receive
#'   new users with complete rankings at each time step. See Chapter 4 of
#'   \insertCite{steinSequentialInferenceMallows2023}{BayesMallows}
#'
#' @param rankings Matrix containing the full set of observed rankings of size
#'   n_assessors by n_items.
#' @param type One of \code{"complete"}, \code{"partial"}, or
#'   \code{"partial_alpha_fixed"}.
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @param leap_size leap_size Integer specifying the step size of the
#'   leap-and-shift proposal distribution
#' @param N Integer specifying the number of particles
#' @param timesteps Integer specifying the number of time steps in the SMC algorithm
#' @param logz_estimate Estimate of the partition function, computed with
#'   \code{\link{estimate_partition_function}}. Be aware that when using an
#'   estimated partition function when \code{n_clusters > 1}, the partition
#'   function should be estimated over the whole range of \eqn{\alpha} values
#'   covered by the prior distribution for \eqn{\alpha} with high probability.
#'   In the case that a cluster \eqn{\alpha_c} becomes empty during the
#'   Metropolis-Hastings algorithm, the posterior of \eqn{\alpha_c} equals its
#'   prior. For example, if the rate parameter of the exponential prior equals,
#'   say \eqn{\lambda = 0.001}, there is about 37 \% (or exactly: \code{1 -
#'   pexp(1000, 0.001)}) prior probability that \eqn{\alpha_c > 1000}. Hence
#'   when \code{n_clusters > 1}, the estimated partition function should cover
#'   this range, or \eqn{\lambda} should be increased.
#' @param mcmc_kernel_app Integer value for the number of applications we apply
#'   the MCMC move kernel
#' @param num_new_obs Integer value for the number of new observations (complete
#'   rankings) for each time step
#' @param alpha_prop_sd Numeric value specifying the standard deviation of the
#'   lognormal proposal distribution used for \eqn{\alpha} in the
#'   Metropolis-Hastings algorithm. Defaults to \code{0.1}.
#' @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the truncated exponential prior distribution of \eqn{\alpha}. Defaults
#'   to \code{0.1}. When \code{n_cluster > 1}, each mixture component
#'   \eqn{\alpha_{c}} has the same prior distribution.
#' @param alpha_max Maximum value of \code{alpha} in the truncated exponential
#'   prior distribution.
#' @param alpha A numeric value of the scale parameter which is known and fixed.
#' @param aug_method A character string specifying the approach for filling in
#'   the missing data, options are "pseudolikelihood" or "random".
#' @param verbose Logical specifying whether to print out the progress of the
#'   SMC-Mallows algorithm. Defaults to \code{FALSE}.
#'
#' @return
#' An object of class \code{"SMCMallows"} containing the following elements:
#' \itemize{
#' \item{"rho_samples"}{An array of samples of the consensus ranking \eqn{\rho}.}
#' \item{"alpha_samples"}{A matrix with samples of \eqn{\alpha}. Empty when \code{alpha_fixed = TRUE}.}
#' \item{"augmented_rankings"}{Array of augmented rankings.}
#' \item{"ESS"}{Vector of effective sample sizes at each timepoint.}
#' }
#'
#' @export
#'
#' @example inst/examples/smc_mallows_new_users_complete_example.R
#'
#' @family modeling
#'
#'
smc_mallows_new_users <- function(
    rankings,
    type = c("complete", "partial", "partial_alpha_fixed"),
    N,
    timesteps,
    mcmc_kernel_app,
    num_new_obs,
    alpha_prop_sd = 0.5,
    lambda = 0.1,
    alpha_max = 1e6,
    alpha = 0,
    aug_method = "random",
    logz_estimate = NULL,
    verbose = FALSE,
    metric = "footrule",
    leap_size = 1) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley", "hamming",
    "kendall", "ulam"
  ))
  stopifnot(is.matrix(rankings))
  type <- match.arg(type, c("complete", "partial", "partial_alpha_fixed"))
  n_items <- ncol(rankings)
  logz_list <- prepare_partition_function(logz_estimate, metric, n_items)
  n_users <- nrow(rankings)

  if (type == "complete" && (timesteps > n_users / num_new_obs)) {
    stop("timesteps should not exceed n_users / num_new_obs.")
  }

  smc_mallows_new_users_cpp(
    rankings = rankings,
    type = type,
    n_items = n_items,
    n_users = n_users,
    N,
    timesteps,
    mcmc_kernel_app,
    num_new_obs,
    alpha_prop_sd,
    lambda,
    alpha_max,
    alpha,
    aug_method,
    logz_estimate = logz_list$logz_estimate,
    cardinalities = logz_list$cardinalities,
    verbose,
    metric,
    leap_size
  )
}
