#' Assessing convergence of Metropolis-Hastings algorithm
#'
#' @param R A matrix of ranked items.
#' @param metric The distance metric to use.
#' @param lambda Parameter for the prior distribution of \code{alpha}. Defaults
#'   to 0.1.
#' @param nmc Number of Monte Carlo samples in test run.
#' @param burnin How many should we regard as burnin? This one is relevant for
#'   computing the acceptance rate.
#' @param L Step size of the leap-and-shift proposal distribution.
#' @param sd_alpha Standard deviation of the proposal distribution for alpha.
#' @param alpha_init Initial value of alpha.
#' @param max_items The maximum number of items to plot in the diagnostic for
#'   rho.
#' @param rho_displays The maximum number of displays in the diagnostic plots
#'   for rho.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @return A list containing two plots, one for alpha and one for rho.
#' @details After finding reasonable values of the tuning parameters, you
#'   typically want to use \code{\link{compute_mallows}}.
#' @export
#'
#' @examples
#' # Let us first look at the convergence of alpha, for the potato data
#' # First we try with an alpha standard deviation of 5
#' check <- assess_convergence(R = potato_weighing, sd_alpha = 5)
#' # Now show the convergence plot for alpha
#' check$alpha_plot
#' # The acceptance rate was a bit low, so let us try to decrease the sd_alpha
#' check <- assess_convergence(R = potato_weighing, sd_alpha = 0.1)
#' # Then view the plot
#' check$alpha_plot
#' # This looked better!
#'
#' # Next, we assess the convergence for rho. The default is L = n/5,
#' # which equals 4 in this case.
#' # We start by using the results of the previous run.
#' check$rho_plot
#' # Let us see how increasing L changes the convergence
#' check <- assess_convergence(R = potato_weighing, sd_alpha = 0.1, L = 1)
#' # Plot the new result
#' check$rho_plot
#' @import rlang
assess_convergence <- function(R, metric = "footrule", lambda = 0.1,
                               nmc = 3000, burnin = 2000,
                               L = ncol(R) / 5, sd_alpha = 0.1, alpha_init = 1,
                               max_items = 20, rho_displays = 4){

  model_fit <- compute_mallows(R, metric, lambda, nmc, burnin = 0, L, sd_alpha, alpha_init)

  # Converting to data frame before plotting
  df <- data.frame(
    index = seq_along(model_fit$alpha),
    alpha = model_fit$alpha,
    alpha_acceptance = model_fit$alpha_acceptance
  )

  alpha_acceptance_rate <- mean(model_fit$alpha_acceptance[(burnin + 1) : nmc])

  # Create the diagnostic plot for alpha
  alpha_plot <- ggplot2::ggplot(df, ggplot2::aes_(x =~ index, y =~ alpha)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = burnin, linetype = "dashed") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha)) +
    ggplot2::ggtitle(
      paste(
        "Acceptance rate after burnin:",
        sprintf("%.1f", alpha_acceptance_rate * 100), "%"
        )
      )

  # Create the diagnostic plot for rho
  # First create a tidy data frame for holding the data
  df <- data.frame(
    Item = as.integer(seq(from = 1, to = nrow(model_fit$rho), by = 1)),
    Iteration = as.integer(rep(seq(from = 1, to = nmc, by = 1), each = nrow(model_fit$rho))),
    Rank = as.integer(model_fit$rho),
    Acceptance = rep(as.integer(model_fit$rho_acceptance), each = nrow(model_fit$rho))
    )

  if(max(df$Item) > max_items) {
    message(paste("More than", max_items, "items. Selecting", max_items, " at random."))
    df <- df[df$Item %in% sample(max(df$Item), max_items), ]
  }

  # We want to plot no more than 5 items per display, so creating an indicator variable
  # for this
  df$Display <- ggplot2::cut_number(df$Item, n = rho_displays,
                                    labels = seq(from = 1, to = rho_displays))
  df$Item <- as.factor(df$Item)

  rho_acceptance <- mean(model_fit$rho_acceptance[(burnin + 1) : nmc])

  rho_plot <- ggplot2::ggplot(df, ggplot2::aes_(x =~ Iteration, y =~ Rank, color =~ Item)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = burnin, linetype = "dashed") +
    ggplot2::facet_wrap(~ Display) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(rho)) +
    ggplot2::ggtitle(
      paste(
        "Acceptance rate after burnin:",
        sprintf("%.1f", rho_acceptance * 100), "%"
      )
    )

  return(list(alpha_plot = alpha_plot, rho_plot = rho_plot))

}
