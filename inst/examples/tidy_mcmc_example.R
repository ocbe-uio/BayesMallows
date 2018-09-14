# Run compute_mallows with the option skip_postprocessing=TRUE.
n_clusters <- 1:10
library(purrr)
models <- n_clusters %>%
  map(~ compute_mallows(potato_visual, n_clusters = .x,
                        include_wcd = TRUE, skip_postprocessing = TRUE))

# We are now not able to compute an elbow plot
\dontrun{
  plot_elbow(models, burnin = 500)
}

# We need to tidy the parts that are used by plot_elbow.
tidy_models <- models %>%
  map(~ tidy_mcmc(.x,
                  tidy_rho = FALSE,
                  tidy_alpha = FALSE,
                  tidy_cluster_assignment = FALSE,
                  tidy_cluster_probabilities = FALSE,
                  tidy_wcd = TRUE,
                  tidy_augmented_data = TRUE,
                  tidy_augmentation_acceptance = TRUE)
      )
plot_elbow(tidy_models, burnin = 500)
