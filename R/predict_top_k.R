#' Predict Top-k Rankings with Pairwise Preferences
#'
#' Predict the posterior probability, per item, of being ranked among the
#' top-\eqn{k} for each assessor. This is useful when the data take the form of
#' pairwise preferences.
#'
#' @param model_fit An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param k Integer specifying the k in top-\eqn{k}.
#'
#' @export
#'
#' @return A dataframe with columns `assessor`, `item`, and
#'   `prob`, where each row states the probability that the given assessor
#'   rates the given item among top-\eqn{k}.
#'
#' @example /inst/examples/plot_top_k_example.R
#' @family posterior quantities
predict_top_k <- function(model_fit,k = 3) {
  validate_top_k(model_fit)
  .predict_top_k(model_fit, k)
}

.predict_top_k <- function(model_fit, k) {
  rankings <- model_fit$augmented_data[
    model_fit$augmented_data$iteration > burnin(model_fit) &
      model_fit$augmented_data$value <= k, ,
    drop = FALSE
  ]

  n_samples <- length(unique(rankings$iteration))
  rankings <- aggregate(
    list(prob = rankings$iteration),
    by = list(assessor = rankings$assessor, item = rankings$item),
    FUN = function(x) length(x) / n_samples, drop = FALSE
  )
  rankings$prob[is.na(rankings$prob)] <- 0
  rankings
}


validate_top_k <- function(model_fit) {
  if (is.null(burnin(model_fit))) {
    stop("Please specify the burnin with 'burnin(model_fit) <- value'.")
  }

  if (!exists("augmented_data", model_fit)) {
    stop("model_fit must have element augmented_data. Please set save_aug = TRUE
         in compute_mallows in order to create a top-k plot.")
  }
}
