#' @title SMC-Mallows new item rank
#' @description Function to perform resample-move SMC algorithm where we receive
#'   a new item ranks from an existing user at each time step. Each correction
#'   and augmentation is done by filling in the missing item ranks using
#'   pseudolikelihood augmentation.
#' @param rankings 3D array of size n_assessors by n_items by timesteps
#'   containing a set of observed rankings of timesteps time steps
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @param leap_size leap_size Integer specifying the step size of the
#'   leap-and-shift proposal distribution
#' @param n_particles Integer specifying the number of particles
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
#' @param mcmc_steps Integer value for the number of applications we apply
#'   the MCMC move kernel
#' @param alpha_prop_sd Numeric value of the standard deviation of the prior
#'   distribution for alpha
#' @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the truncated exponential prior distribution of alpha.
#' @param alpha_max  Maximum value of alpha in the truncated exponential prior
#'   distribution.
#' @param aug_method A character string specifying the approach for filling in
#'   the missing data, options are "pseudolikelihood" or "random"
#' @param verbose Logical specifying whether to print out the progress of the
#'   SMC-Mallows algorithm. Defaults to \code{FALSE}.
#' @param alpha_fixed Logical indicating whether to sample \code{alpha} or not.
#' @param alpha numeric value of the scale parameter.
#' @param aug_rankings_init Initial values for augmented rankings.
#' @param rho_samples_init Initial values for rho samples.
#' @param alpha_samples_init Initial values for alpha samples.
#'
#' @return An object of class \code{c("SMCMallowsUpdatedPartial", "SMCMallows")}
#' containing the following elements:
#' \itemize{
#' \item{"rho_samples"}{An array of samples of the consensus ranking \eqn{\rho}.}
#' \item{"alpha_samples"}{A matrix with samples of \eqn{\alpha}. Empty when \code{alpha_fixed = TRUE}.}
#' \item{"augmented_rankings"}{Array of augmented rankings.}
#' \item{"ESS"}{Vector of effective sample sizes at each timepoint.}
#' }
#'
#' @export
#'
#' @family modeling
#'
#'
smc_mallows_new_item_rank <- function(
    rankings,
    n_particles,
    logz_estimate = NULL,
    mcmc_steps,
    aug_rankings_init = NULL,
    rho_samples_init = NULL,
    alpha_samples_init = 0,
    alpha = 0,
    alpha_prop_sd = 0.5,
    lambda = 0.1,
    alpha_max = 1e6,
    aug_method = "random",
    verbose = FALSE,
    alpha_fixed = FALSE,
    metric = "footrule",
    leap_size = 1) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley", "hamming",
    "kendall", "ulam"
  ))

  n_items <- dim(rankings)[[2]]
  timesteps <- dim(rankings)[[3]]

  logz_list <- prepare_partition_function(logz_estimate, metric, n_items)

  ret <- smc_mallows_new_item_rank_cpp(
    n_items,
    rankings,
    n_particles,
    timesteps,
    logz_estimate = logz_list$logz_estimate,
    cardinalities = logz_list$cardinalities,
    mcmc_steps,
    aug_rankings_init,
    rho_samples_init,
    alpha_samples_init,
    alpha,
    alpha_prop_sd,
    lambda,
    alpha_max,
    aug_method,
    verbose,
    alpha_fixed,
    metric,
    leap_size
  )
  ret$n_particles <- n_particles
  ret$metric <- metric
  ret$logz_list <- logz_list
  ret$n_items <- n_items
  ret$mcmc_steps <- mcmc_steps
  ret$alpha_prop_sd <- alpha_prop_sd
  ret$lambda <- lambda
  ret$alpha_max <- alpha_max
  ret$alpha <- alpha
  ret$alpha_fixed <- alpha_fixed
  ret$aug_method <- aug_method
  ret$leap_size <- leap_size
  ret$num_obs <- dim(rankings)[[1]]
  ret$rankings <- rankings
  class(ret) <- c("SMCMallowsUpdatedPartial", "SMCMallows")
  ret
}
