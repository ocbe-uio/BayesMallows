#' @title Set the burnin
#' @description Set or update the burnin of a model
#'   computed using Metropolis-Hastings.
#'
#' @param model An object of class `BayesMallows` returned from
#'   [compute_mallows()] or an object of class `BayesMallowsMixtures` returned
#'   from [compute_mallows_mixtures()].
#' @param ... Optional arguments passed on to other methods. Currently not used.
#' @param value An integer specifying the burnin. If `model` is of class
#'   `BayesMallowsMixtures`, a single value will be assumed to be the burnin
#'   for each model element. Alternatively, `value` can be specified as an
#'   integer vector of the same length as `model`, and hence a separate burnin
#'   can be set for each number of mixture components.
#'
#' @export
#' @return An object of class `BayesMallows` with burnin set.
#'
#' @family modeling
#'
#' @example /inst/examples/burnin_example.R
#'
`burnin<-` <- function(model, ..., value) UseMethod("burnin<-")

#' @export
#' @rdname burnin-set
`burnin<-.BayesMallows` <- function(model, ..., value) {
  if (inherits(model, "SMCMallows")) {
    stop("Cannot set burnin for SMC model.")
  }
  validate_integer(value)
  if (value >= model$compute_options$nmc) {
    stop("Burnin cannot be larger than the number of Monte Carlo samples.")
  }
  # Workaround as long as we have the deprecation notice for `$<-`
  class(model) <- "list"
  model$compute_options$burnin <- value
  class(model) <- "BayesMallows"
  model
}

#' @export
#' @rdname burnin-set
`burnin<-.BayesMallowsMixtures` <- function(model, ..., value) {
  for (v in value) validate_integer(v)
  if (length(value) == 1) value <- rep(value, length(model))
  if (length(value) != length(model)) stop("Wrong number of entries in value.")

  for (i in seq_along(model)) burnin(model[[i]]) <- value[[i]]
  model
}

#' @title See the burnin
#' @description
#' See the current burnin value of the model.
#'
#' @param model A model object.
#' @param ... Optional arguments passed on to other methods. Currently not used.
#'
#' @export
#' @return An integer specifying the burnin, if it exists. Otherwise `NULL`.
#'
#' @family modeling
#'
#' @example /inst/examples/burnin_example.R
#'
burnin <- function(model, ...) UseMethod("burnin")

#' @rdname burnin
#' @export
burnin.BayesMallows <- function(model, ...) {
  model$compute_options$burnin
}

#' @rdname burnin
#' @export
burnin.BayesMallowsMixtures <- function(model, ...) {
  lapply(model, burnin)
}

#' @rdname burnin
#' @export
burnin.SMCMallows <- function(model, ...) 0
