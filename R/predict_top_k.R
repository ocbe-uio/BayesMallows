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
predict_top_k <- function(model_fit, burnin = model_fit$burnin,
                          k = 3){

  validate_top_k(model_fit, burnin)
  .predict_top_k(model_fit, burnin, k)
}



.predict_top_k <- function(model_fit, burnin, k){

  rankings <- dplyr::filter(model_fit$augmented_data, .data$iteration > burnin, .data$value <= k)
  n_samples <- length(unique(rankings$iteration))
  rankings <- dplyr::mutate(rankings, item = as.character(.data$item))
  rankings <- dplyr::group_by(rankings, .data$assessor, .data$item)
  rankings <- dplyr::summarise(rankings, prob = dplyr::n()/n_samples, .groups = "drop")

  rankings <- tidyr::complete(
    dplyr::group_by(rankings, .data$assessor),
    item = model_fit$items,
    fill = list(prob = 0)
  )


  return(rankings)
}


validate_top_k <- function(model_fit, burnin){
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  if(!exists("augmented_data", model_fit)){
    stop("model_fit must have element augmented_data. Please set save_aug = TRUE
         in compute_mallows in order to create a top-k plot.")
  }
  }
