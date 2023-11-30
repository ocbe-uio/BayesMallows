library(tidyverse)
# We use the convenience function compute_mallows_mixtures
n_clusters <- c(4, 8)
models <- compute_mallows_mixtures(
  n_clusters = n_clusters, data = setup_rank_data(sushi_rankings))
# models is a list in which each element is an object of class BayesMallows,
# returned from compute_mallows
# We can create an elbow plot
plot_elbow(models, burnin = 1000)
# We then select the number of cluster at a point where this plot has
# an "elbow", e.g., n_clusters = 5.

# Having chosen the number of clusters, we can now study the final model
# Rerun with 5 clusters
mixture_model <- compute_mallows(
  data = setup_rank_data(rankings = sushi_rankings),
  compute_options = set_compute_options(nmc = 700),
  model_options = set_model_options(n_clusters = 5))


apply(mixture_model$alpha_samples, 1, mean)

mixture_model$cluster_assignment %>%
  as_tibble() %>%
  group_by(value) %>%
  summarise(n())

compute_consensus(mixture_model, burnin = 200)
plot(mixture_model, parameter = "cluster_probs", burnin = 200)
