#' Set options for Bayesian Mallows model
#'
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}. The distance given by \code{metric} is also used to compute
#'   within-cluster distances, when \code{include_wcd = TRUE}.
#'
#'
#' @param error_model Character string specifying which model to use for
#'   inconsistent rankings. Defaults to \code{"none"}, which means that
#'   inconsistent rankings are not allowed. At the moment, the only available
#'   other option is \code{"bernoulli"}, which means that the Bernoulli error
#'   model is used. See \insertCite{crispino2019;textual}{BayesMallows} for a
#'   definition of the Bernoulli model.
#'
#'
#' @param n_clusters Integer specifying the number of clusters, i.e., the number
#'   of mixture components to use. Defaults to \code{1L}, which means no
#'   clustering is performed. See \code{\link{compute_mallows_mixtures}} for a
#'   convenience function for computing several models with varying numbers of
#'   mixtures.
#'
#'
#' @return An object of class \code{"BayesMallowsModelOptions"}, to be provided
#'   in the \code{model} argument to \code{\link{compute_mallows}}.
#'
#' @export
#'
#' @family options
#'
set_model_options <- function(metric = "footrule", n_clusters = 1,
                              error_model = "none") {
  metric <- match.arg(metric, c("footrule", "spearman", "cayley", "hamming",
                                "kendall", "ulam"))
  error_model <- match.arg(error_model, c("none", "bernoulli"))

  validate_integer(n_clusters)
  validate_positive(n_clusters)

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsModelOptions"
  ret
}
