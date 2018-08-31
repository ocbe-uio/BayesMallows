\dontrun{
  # IMPORTANCE SAMPLING
  # Let us estimate logZ(alpha) for 20 items with Spearman distance
  # We create a grid of alpha values from 0 to 10
  alpha_vector <- seq(from = 0, to = 10, by = 0.5)
  n_items <- 20
  metric <- "spearman"
  degree <- 10

  # We start with 1e3 Monte Carlo samples
  fit1 <- estimate_partition_function(method = "importance_sampling",
                                        alpha_vector = alpha_vector,
                                        n_items = n_items, metric = metric,
                                        nmc = 1e3, degree = degree)
  # A vector of polynomial regression coefficients is returned
  fit1

  # Now let us recompute with 1e4 Monte Carlo samples
  fit2 <- estimate_partition_function(method = "importance_sampling",
                                      alpha_vector = alpha_vector,
                                      n_items = n_items, metric = metric,
                                      nmc = 1e4, degree = degree)

  # ASYMPTOTIC APPROXIMATION
  # We can also compute an estimate using the asymptotic approximation
  K <- 20
  n_iterations <- 50

  fit3 <- estimate_partition_function(method = "asymptotic",
                                      alpha_vector = alpha_vector,
                                      n_items = n_items, metric = metric,
                                      n_iterations = n_iterations,
                                      K = K, degree = degree)

  # We write a little function for storing the estimates in a dataframe
  library(dplyr)
  powers <- seq(from = 0, to = degree, by = 1)

  compute_fit <- function(fit){
    tibble(alpha = alpha_vector) %>%
      rowwise() %>%
      mutate(logz_estimate = sum(alpha^powers * fit))
  }

  estimates <- bind_rows(
    "Importance Sampling 1e3" = compute_fit(fit1),
    "Importance Sampling 1e4" = compute_fit(fit2),
    "Asymptotic" = compute_fit(fit3),
    .id = "type")

  # We can now plot the two estimates side-by-side
  library(ggplot2)
  ggplot(estimates, aes(x = alpha, y = logz_estimate, color = type)) +
    geom_line()
  # We see that the two importance sampling estimates, which are unbiased,
  # overlap. The asymptotic approximation seems a bit off. It can be worthwhile
  # to try different values of n_iterations and K.

  # When we are happy, we can provide the coefficient vector in the
  # logz_estimate argument to compute_mallows
  # Say we choose to use the importance sampling estimate with 1e4 Monte Carlo samples:
  model_fit <- compute_mallows(potato_visual, metric = "spearman",
                               logz_estimate = fit2)

}
