# DETERMINING THE NUMBER OF CLUSTERS IN THE SUSHI EXAMPLE DATA
\dontrun{
  # Let us look at any number of clusters from 1 to 10
  # We use the convenience function compute_mallows_mixtures
  n_clusters <- seq(from = 1, to = 10)
  models <- compute_mallows_mixtures(n_clusters = n_clusters,
                                     rankings = sushi_rankings,
                                     include_wcd = TRUE)
  # models is a list in which each element is an object of class BayesMallows,
  # returned from compute_mallows
  # We can create an elbow plot
  plot_elbow(models, burnin = 1000)
  # We then select the number of cluster at a point where this plot has
  # an "elbow", e.g., n_clusters = 5.

  # Having chosen the number of clusters, we can now study the final model
  # Rerun with 5 clusters, now setting save_clus = TRUE to get cluster assignments
  mixture_model <- compute_mallows(rankings = sushi_rankings, n_clusters = 5,
                                   include_wcd = TRUE, save_clus = TRUE)
  # Delete the models object to free some memory
  rm(models)
  # Set the burnin
  mixture_model$burnin <- 1000
  # Plot the posterior distributions of alpha per cluster
  plot(mixture_model)
  # Compute the posterior interval of alpha per cluster
  compute_posterior_intervals(mixture_model,
                              parameter = "alpha")
  # Plot the posterior distributions of cluster probabilities
  plot(mixture_model, parameter = "cluster_probs")
  # Plot the posterior probability of cluster assignment
  plot(mixture_model, parameter = "cluster_assignment")
  # Plot the posterior distribution of "tuna roll" in each cluster
  plot(mixture_model, parameter = "rho", items = "tuna roll")
  # Compute the cluster-wise CP consensus, and show one column per cluster
  cp <- compute_consensus(mixture_model, type = "CP")
  library(dplyr)
  library(tidyr)
  cp %>%
    select(-cumprob) %>%
    spread(key = cluster, value = item)
  # Compute the MAP consensus, and show one column per cluster
  map <- compute_consensus(mixture_model, type = "MAP")
  map %>%
    select(-probability) %>%
    spread(key = cluster, value = item)

  # RUNNING IN PARALLEL
  # Computing Mallows models with different number of mixtures in parallel leads to
  # considerably speedup
  library(parallel)
  cl <- makeCluster(detectCores() - 1)
  n_clusters <- seq(from = 1, to = 10)
  models <- compute_mallows_mixtures(n_clusters = n_clusters,
                                     rankings = sushi_rankings,
                                     include_wcd = TRUE, cl = cl)
  stopCluster(cl)
}



