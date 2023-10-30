#' Compute Mixtures of Mallows Models
#'
#' Convenience function for computing Mallows models with varying
#' numbers of mixtures. This is useful for deciding the number of
#' mixtures to use in the final model.
#'
#' @param n_clusters Integer vector specifying the number of clusters
#' to use.
#'
#' @param ... Other named arguments, passed to \code{\link{compute_mallows}}.
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
#'
#' @return A list of Mallows models of class \code{BayesMallowsMixtures}, with one element
#' for each number of mixtures that
#' was computed. This object can be studied with \code{\link{plot_elbow}}.
#'
#' @family modeling
#'
#' @export
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#'
compute_mallows_mixtures <- function(n_clusters, ..., cl = NULL) {
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  call <- match.call()
  if(is.null(call$model)) call$model <- set_model_options()
  call[[1]] <- substitute(compute_mallows)
  call[["n_clusters"]] <- NULL

  if (is.null(cl)) {
    lapplyfun <- lapply
  } else {

    find_objects_to_export <- function(arg) {
      if(is.call(arg)) {
        lapply(arg[-1], find_objects_to_export)
      } else if(is.name(arg)) {
        arg
      } else {
        character()
      }
    }

    call[["cl"]] <- NULL
    args <- as.list(call[-1])

    varlist <- as.character(unlist(lapply(args, find_objects_to_export)))
    varlist <- varlist[varlist != ""]
    parallel::clusterExport(cl = cl, varlist = c("call", varlist),
                            envir = environment())
    lapplyfun <- function(X, FUN, ...) {
      parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    }
  }

  models <- lapplyfun(n_clusters, function(x) {
    call$model$n_clusters <- x
    eval(call)
  })

  class(models) <- "BayesMallowsMixtures"
  return(models)
}
