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

    df <- prepare_alpha_df(model_fit$alpha)

    # Create the diagnostic plot for alpha

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$index, y = .data$alpha)) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(alpha)) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(label = "Convergence of alpha")

    if(model_fit$n_clusters == 1){
      p <- p + ggplot2::geom_line()
    } else {
      p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$cluster))
    }

    print(p)

  } else if(type == "rho"){

    if(is.null(items) && model_fit$n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(model_fit$n_items, 5)
    } else if (is.null(items) && model_fit$n_items > 0) {
      items <- seq.int(from = 1, to = model_fit$n_items)
    }

    df <- gather_rho(model_fit, items)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$index, y = .data$rank, color = .data$item)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(rho)) +
      ggplot2::ggtitle(label = "Convergence of rho")

    if(model_fit$n_clusters > 1){
      p <- p + ggplot2::facet_wrap(~ .data$cluster)
    }

    print(p)

  } else if(type == "augmentation") {
    if(!model_fit$any_missing && !model_fit$augpair) stop("No missing values, so data were not augmented.")

    if(is.null(assessors) && model_fit$n_assessors > 12) {
      message("Assessors not set Plotting the first 12.")
      assessors <- seq(from = 1, to = 12, by = 1)
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
                          model_fit$aug_diag_thinning)

    # Saving and then printing, because one gets a warning due to log(0) == -Inf
    ggplot2::ggplot(df, ggplot2::aes(x = .data$Iteration, y = .data$Acceptance)) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::facet_wrap(~ .data$Assessor) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Mean acceptance rate") +
      ggplot2::ggtitle(
        label = "Data Augmentation",
        subtitle = paste("Average per", model_fit$aug_diag_thinning, "steps."))


  } else {
    stop("type must be either \"alpha\", \"rho\", or \"augmentation\".")
  }
}
