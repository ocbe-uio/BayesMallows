\dontrun{
  # We use the example dataset with beach preferences. Se the documentation to
  # compute_mallows for how to assess the convergence of the algorithm
  # We need to save the augmented data, so setting this option to TRUE
  model_fit <- compute_mallows(preferences = beach_preferences,
                               save_aug = TRUE)
  # We set burnin = 1000
  model_fit$burnin <- 1000
  # By default, the probability of being top-3 is plotted
  plot_top_k(model_fit)
  # We can also plot the probability of being top-5, for each item
  plot_top_k(model_fit, k = 5)
  # We get the underlying numbers with predict_top_k
  probs <- predict_top_k(model_fit)
  # To find all items ranked top-3 by assessors 1-3 with probability more than 80 %,
  # we do
  library(dplyr)
  probs %>%
    filter(assessor %in% 1:3, prob > 0.8)

}
