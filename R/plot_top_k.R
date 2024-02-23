#' Plot Top-k Rankings with Pairwise Preferences
#'
#' Plot the posterior probability, per item, of being ranked among the
#' top-\eqn{k} for each assessor. This plot is useful when the data take the
#' form of pairwise preferences.
#'
#' @param model_fit An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param k Integer specifying the k in top-\eqn{k}.
#'
#' @export
#'
#' @example /inst/examples/plot_top_k_example.R
#' @family posterior quantities
plot_top_k <- function(model_fit, k = 3) {
  validate_top_k(model_fit)
  rankings <- .predict_top_k(model_fit, k = k)
  ggplot2::ggplot(rankings, ggplot2::aes(.data$assessor, .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "white") +
    ggplot2::xlab("Assessor") +
    ggplot2::ylab("Item") +
    ggplot2::labs(fill = "Prob.")
}
