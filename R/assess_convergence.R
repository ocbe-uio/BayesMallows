#' Create plots for assessing convergence of Metropolis-Hastings algorithm
#'
#' @param R Ranks
#' @param metric string
#' @param lambda Prior
#' @param nmc Number of Monte Carlo samples in test run
#' @param burnin How many should we regard as burnin?
#' @param L Leap-and-shift step size
#' @param sd_alpha Standard deviation of alpha proposal
#' @param alpha_init Initial value of alpha
#'
#' @return A list containing two plots, one for alpha and one for rho.
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
assess_convergence <- function(R, metric = "footrule", lambda = 0.1,
                               nmc = 10000, burnin = 5000,
                               L = 1, sd_alpha = 0.1, alpha_init = 1){
  model_fit <- compute_posterior(R, metric, lambda, nmc, L, sd_alpha, alpha_init)


  alpha_acceptance_rate <- sum(model_fit$alpha_acceptance[(burnin + 1) : nmc]) / (nmc - burnin)

  # Converting to data frame before plotting
  df <- data.frame(
    index = seq_along(model_fit$alpha),
    alpha = model_fit$alpha,
    alpha_acceptance = model_fit$alpha_acceptance
  )

  alpha_plot <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = alpha)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = eval(burnin), linetype = "dashed") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha)) +
    ggplot2::ggtitle(
      paste(
        "Acceptance rate after burnin:",
        sprintf("%.1f", alpha_acceptance_rate * 100), "%"
        )
      )


  return(list(alpha_plot = alpha_plot))

}
