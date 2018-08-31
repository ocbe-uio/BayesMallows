# DETERMINING THE NUMBER OF CLUSTERS
\dontrun{
  # Please see the documentation to compute_mallows for how to assess
  # convergence and plot posterior distributions of parameters.
  # We use the example dataset with sushi preferences.
  # Let us look at any number of clusters from 1 to 10
  # We use the map function from the purrr package to call compute_mallows
  # repeatedly and storing the results in a list named models.
  library(purrr)
  n_clusters <- seq(from = 1, to = 10)
  # Setting the include_wcd argument to TRUE ensures that the within-cluster distances
  # (WCD) are saved in the MCMC algorithm. It defaults to TRUE when n_clusters > 1,
  # but we also need it for n_clusters = 1, and hence set it explicitly.
  models <- map(n_clusters, ~ compute_mallows(rankings = sushi_rankings, nmc = 6000,
                                              n_clusters = .x, include_wcd = TRUE))
  # models is a list in which each element is an object of class BayesMallows,
  # returned from compute_mallows
  # We can create an elbow plot
  plot_elbow(models, burnin = 1000)
  # We then select the number of cluster at a point where this plot has
  # an "elbow", e.g., n_clusters = 5.
}



