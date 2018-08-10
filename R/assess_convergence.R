#' Assessing convergence of Metropolis-Hastings algorithm
#'
#' @param model_fit A fitted model object of class \code{BayesMallows}, obtained
#'   with \code{\link{compute_mallows}}.
#' @param type Should we plot \code{"alpha"}, \code{"rho"}, or
#'   \code{"augmentation"}.
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{model_fit$item_names} or a
#'   vector of indices. If NULL, five items are selected randomly.
#' @param assessors The assessors to study in the diagnostic plot for
#'   \code{"augmentation"}.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @export
#'
assess_convergence <- function(model_fit, type = "alpha", items = NULL,
                               assessors = NULL){

  stopifnot(class(model_fit) == "BayesMallows")

  if(is.character(items)) stopifnot(items %in% rownames(model_fit$rho))

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

  } else if(type == "augmentation") {
    if(!model_fit$any_missing) stop("No missing values, so data were not augmented.")

    if(is.null(assessors) && model_fit$n_assessors > 12) {
      message("Assessors not set Plotting the first 12.")
    } else if(is.null(assessors)) {
      assessors <- seq(from = 1, to = model_fit$n_assessors, by = 1)
    }

    df <- dplyr::as_tibble(model_fit$aug_acceptance)
    names(df) <- 1:ncol(df)

    df <- dplyr::mutate(df, Assessor = as.factor(dplyr::row_number()))
    df <- dplyr::slice(df, assessors)

    df <- tidyr::gather(df, key = "Iteration", value = "Acceptance",
                        -.data$Assessor, convert = TRUE)

    df <- dplyr::mutate(df, Iteration = .data$Iteration *
                          model_fit$augmentation_diagnostic_thinning)

    # Saving and then printing, because one gets a warning due to log(0) == -Inf
    ggplot2::ggplot(df, ggplot2::aes(x = .data$Iteration, y = .data$Acceptance)) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::facet_wrap(~ .data$Assessor) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Mean acceptance rate") +
      ggplot2::ggtitle(
        label = "Data Augmentation",
        subtitle = paste("Average per", model_fit$augmentation_diagnostic_thinning, "steps."))


  } else {
    stop("type must be either \"alpha\", \"rho\", or \"augmentation\".")
  }
}
