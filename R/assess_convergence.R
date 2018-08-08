#' Assessing convergence of Metropolis-Hastings algorithm
#'
#' @param model_fit A fitted model object of class \code{BayesMallows}, obtained
#'   with \code{\link{compute_mallows}}.
#' @param items The items to study in the diagnostic plot for \code{rho}. A
#'   vector of indices. If NULL, five items are selected randomly.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @return A list containing two plots, one for alpha and one for rho.
#' @export
#'
assess_convergence <- function(model_fit, items = NULL){

  stopifnot(class(model_fit) == "BayesMallows")

  # Converting to data frame before plotting
  df <- dplyr::bind_cols(
    alpha = model_fit$alpha,
    alpha_acceptance = model_fit$alpha_acceptance)

  df <- dplyr::mutate(df, Index = dplyr::row_number())

  alpha_acceptance_rate <- dplyr::pull(dplyr::summarise(df, mean(.data$alpha_acceptance)))

  # Create the diagnostic plot for alpha
  alpha_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$alpha)) +
    ggplot2::geom_line() +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha)) +
    ggplot2::ggtitle(
      label = "Convergence of alpha",
      subtitle = paste("Acceptance rate:", sprintf("%.1f", alpha_acceptance_rate * 100), "%")
      )

  # Create the diagnostic plot for rho
  # First create a tidy data frame for holding the data

  if(is.null(items) && model_fit$n_items > 5){
    items <- sample.int(model_fit$n_items, 5)
  } else if (is.null(items) && model_fit$n_items > 0) {
    items <- seq.int(from = 1, to = model_fit$n_items)
  }

  df <- gather_rho(model_fit, items)

  rho_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$Rank, color = .data$Item)) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(rho)) +
    ggplot2::ggtitle(label = "Convergence of rho")

  return(list(alpha_plot = alpha_plot, rho_plot = rho_plot))

}
