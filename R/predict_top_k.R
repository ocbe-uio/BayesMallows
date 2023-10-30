#' Predict Top-k Rankings with Pairwise Preferences
#'
#' Predict the posterior probability, per item, of being ranked among the
#' top-\eqn{k} for each assessor. This is useful when the data take the form of
#' pairwise preferences.
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations to discard
#'   as burn-in. Defaults to \code{model_fit$burnin}, and must be provided if
#'   \code{model_fit$burnin} does not exist. See
#'   \code{\link{assess_convergence}}.
#'
#' @param k Integer specifying the k in top-\eqn{k}.
#'
#' @export
#'
#' @return A dataframe with columns \code{assessor}, \code{item}, and
#'   \code{prob}, where each row states the probability that the given assessor
#'   rates the given item among top-\eqn{k}.
#'
#' @seealso \code{\link{plot_top_k}}
#'
#' @example /inst/examples/plot_top_k_example.R
#' @family posterior quantities
predict_top_k <- function(model_fit, burnin = model_fit$burnin, k = 3) {
    validate_top_k(model_fit, burnin)
    .predict_top_k(model_fit, burnin, k)
}

.predict_top_k <- function(model_fit, burnin, k) {
  rankings <- model_fit$augmented_data[model_fit$augmented_data$iteration > burnin &
    model_fit$augmented_data$value <= k, , drop = FALSE]

  n_samples <- length(unique(rankings$iteration))
  rankings$item <- as.character(rankings$item)
  rankings <- aggregate(
    list(prob = rankings$iteration),
    by = list(assessor = rankings$assessor, item = rankings$item),
    FUN = function(x) length(x) / n_samples, drop = FALSE
  )
  rankings$prob[is.na(rankings$prob)] <- 0

  rankings[order(rankings$assessor, rankings$item), ]
}


validate_top_k <- function(model_fit, burnin) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  if (!exists("augmented_data", model_fit)) {
    stop("model_fit must have element augmented_data. Please set save_aug = TRUE
         in compute_mallows in order to create a top-k plot.")
  }
}
