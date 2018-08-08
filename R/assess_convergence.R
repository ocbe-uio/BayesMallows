#' Assessing convergence of Metropolis-Hastings algorithm
#'
#' @param model_fit A fitted model object of class \code{BayesMallows}, obtained
#'   with \code{\link{compute_mallows}}.
#' @param type Should we plot \code{alpha} or \code{rho}.
#' @param items The items to study in the diagnostic plot for \code{rho}. A
#'   vector of indices. If NULL, five items are selected randomly.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @return A list containing two plots, one for alpha and one for rho.
#' @export
#'
assess_convergence <- function(model_fit, type = "alpha", items = NULL){

  stopifnot(class(model_fit) == "BayesMallows")

  if(type == "alpha") {
    # Converting to data frame before plotting
    df <- dplyr::bind_cols(
      alpha = model_fit$alpha,
      alpha_acceptance = model_fit$alpha_acceptance)

    df <- dplyr::mutate(df, Index = dplyr::row_number())

    alpha_acceptance_rate <- mean(model_fit$alpha_acceptance)

    # Create the diagnostic plot for alpha
    ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$alpha)) +
      ggplot2::geom_line() +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(alpha)) +
      ggplot2::ggtitle(
        label = "Convergence of alpha",
        subtitle = paste("Acceptance rate:",
                         sprintf("%.1f", alpha_acceptance_rate * 100), "%")
      )

  } else if(type == "rho"){
    # Create the diagnostic plot for rho
    # First create a tidy data frame for holding the data

    if(is.null(items) && model_fit$n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(model_fit$n_items, 5)
    } else if (is.null(items) && model_fit$n_items > 0) {
      items <- seq.int(from = 1, to = model_fit$n_items)
    }

    df <- gather_rho(model_fit, items)
    rho_acceptance_rate <- mean(model_fit$rho_acceptance)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$Rank, color = .data$Item)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(rho)) +
      ggplot2::ggtitle(
        label = "Convergence of rho",
        subtitle = paste("Acceptance rate:",
                         sprintf("%.1f", rho_acceptance_rate * 100), "%"))
  } else {
    stop("type must be either \"alpha\" or \"rho\"")
  }
}
