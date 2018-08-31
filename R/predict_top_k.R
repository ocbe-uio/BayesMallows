#' Predict Top-k Rankings with Pairwise Preferences
#'
#' Predict the posterior probability, per item, of being ranked among the top-\eqn{k}
#' for each assessor. This is useful when the data take the form of pairwise
#' preferences. This is an internal function used by \code{\link{plot_top_k}}.
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. See \code{\link{assess_convergence}}.
#'
#' @param k Integer specifying the k in top-\eqn{k}.
#'
#' @param n_samples The number of samples, after burn-in.
#'
#' @keywords internal
#'
#' @seealso \code{\link{plot_top_k}}
#'
predict_top_k <- function(model_fit, burnin, k, n_samples){

  rankings <- dplyr::filter(model_fit$augmented_data, .data$iteration > burnin, .data$value <= k)
  rankings <- dplyr::mutate(rankings, item = as.character(.data$item))
  rankings <- dplyr::group_by(rankings, .data$assessor, .data$item)
  rankings <- dplyr::summarise(rankings, prob = dplyr::n()/n_samples)
  rankings <- dplyr::ungroup(rankings)
  rankings <- tidyr::complete(
    dplyr::group_by(rankings, .data$assessor),
    item = model_fit$items,
    fill = list(prob = 0)
  )


  return(rankings)
}
