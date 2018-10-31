# ANALYSIS OF COMPLETE RANKINGS
# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
# model:
model_fit <- compute_mallows(potato_visual)

# We study the trace plot of the parameters
assess_convergence(model_fit, parameter = "alpha")
\dontrun{assess_convergence(model_fit, parameter = "rho")}

# Based on these plots, we set burnin = 1000.
model_fit$burnin <- 1000
# Next, we use the generic plot function to study the posterior distributions
# of alpha and rho
plot(model_fit, parameter = "alpha")
\dontrun{plot(model_fit, parameter = "rho", items = 10:15)}

# We can also compute the CP consensus posterior ranking
compute_consensus(model_fit, type = "CP")

# And we can compute the posterior intervals:
# First we compute the interval for alpha
compute_posterior_intervals(model_fit, parameter = "alpha")
# Then we compute the interval for all the items
\dontrun{compute_posterior_intervals(model_fit, parameter = "rho")}

# ANALYSIS OF PAIRWISE PREFERENCES
\dontrun{
  # The example dataset beach_preferences contains pairwise
  # preferences between beaches stated by 60 assessors. There
  # is a total of 15 beaches in the dataset.
  # In order to use it, we first generate all the orderings
  # implied by the pairwise preferences.
  beach_tc <- generate_transitive_closure(beach_preferences)
  # We also generate an inital rankings
  beach_rankings <- generate_initial_ranking(beach_tc, n_items = 15)
  # We then run the Bayesian Mallows rank model
  # We save the augmented data for diagnostics purposes.
  model_fit <- compute_mallows(rankings = beach_rankings,
                               preferences = beach_tc,
                               save_aug = TRUE,
                               verbose = TRUE)
  # We can assess the convergence of the scale parameter
  assess_convergence(model_fit)
  # We can assess the convergence of latent rankings. Here we
  # show beaches 1-5.
  assess_convergence(model_fit, parameter = "rho", items = 1:5)
  # We can also look at the convergence of the augmented rankings for
  # each assessor.
  assess_convergence(model_fit, parameter = "Rtilde",
                     items = c(2, 4), assessors = c(1, 2))
  # Notice how, for assessor 1, the lines cross each other, while
  # beach 2 consistently has a higher rank value (lower preference) for
  # assessor 2. We can see why by looking at the implied orderings in
  # beach_tc
  library(dplyr)
  beach_tc %>%
    filter(assessor %in% c(1, 2),
           bottom_item %in% c(2, 4) & top_item %in% c(2, 4))
  # Assessor 1 has no implied ordering between beach 2 and beach 4,
  # while assessor 2 has the implied ordering that beach 4 is preferred
  # to beach 2. This is reflected in the trace plots.
}

# CLUSTERING OF ASSESSORS WITH SIMILAR PREFERENCES
\dontrun{
  # The example dataset sushi_rankings contains 5000 complete
  # rankings of 10 types of sushi
  # We start with computing a 3-cluster solution, and save
  # cluster assignments by setting save_clus = TRUE
  model_fit <- compute_mallows(sushi_rankings, n_clusters = 3,
                               nmc = 10000, save_clus = TRUE, verbose = TRUE)
  # We then assess convergence of the scale parameter alpha
  assess_convergence(model_fit)
  # Next, we assess convergence of the cluster probabilities
  assess_convergence(model_fit, parameter = "cluster_probs")
  # Based on this, we set burnin = 1000
  # We now plot the posterior density of the scale parameters alpha in
  # each mixture:
  model_fit$burnin <- 1000
  plot(model_fit, parameter = "alpha")
  # We can also compute the posterior density of the cluster probabilities
  plot(model_fit, parameter = "cluster_probs")
  # We can also plot the posterior cluster assignment. In this case,
  # the assessors are sorted according to their maximum a posteriori cluster estimate.
  plot(model_fit, parameter = "cluster_assignment")
  # We can also assign each assessor to a cluster
  cluster_assignments <- assign_cluster(model_fit, soft = FALSE)
  }

# DETERMINING THE NUMBER OF CLUSTERS
\dontrun{
  # Continuing with the sushi data, we can determine the number of cluster
  # Let us look at any number of clusters from 1 to 10
  # We use the convenience function compute_mallows_mixtures
  n_clusters <- seq(from = 1, to = 10)
  models <- compute_mallows_mixtures(n_clusters = n_clusters, rankings = sushi_rankings,
                                     nmc = 6000, alpha_jump = 10, include_wcd = TRUE)
  # models is a list in which each element is an object of class BayesMallows,
  # returned from compute_mallows
  # We can create an elbow plot
  plot_elbow(models, burnin = 1000)
  # We then select the number of cluster at a point where this plot has
  # an "elbow", e.g., at 6 clusters.
}

