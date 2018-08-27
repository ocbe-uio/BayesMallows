# ANALYSIS OF COMPLETE RANKINGS
# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
# model:
model_fit <- compute_mallows(potato_visual, nmc = 10000)

# We study the trace plot of the parameters
assess_convergence(model_fit, type = "alpha")
\dontrun{assess_convergence(model_fit, type = "rho")}

# Based on these plots, we set burnin = 1000.
# Next, we use the generic plot function to study the posterior distributions
# of alpha and rho
plot(model_fit, burnin = 1000, type = "alpha")
\dontrun{plot(model_fit, burnin = 1000, type = "rho", items = 10:15)}

# We can also compute the CP consensus posterior ranking
compute_cp_consensus(model_fit, burnin = 1000)

# And we can compute the posterior intervals:
# First we compute the interval for alpha
compute_posterior_intervals(model_fit, burnin = 1000, parameter = "alpha")
# Then we compute the interval for all the items
\dontrun{compute_posterior_intervals(model_fit, burnin = 1000, parameter = "rho")}

# CLUSTERING OF ASSESSORS WITH SIMILAR PREFERENCES
\dontrun{
  # The example dataset sushi_rankings contains 5000 complete rankings of 10 types of sushi
  # We start with computing a 3-cluster solution
  model_fit <- compute_mallows(sushi_rankings, n_clusters = 3, nmc = 10000)
  # We then assess convergence of the scale parameter alpha
  assess_convergence(model_fit)
  # Next, we assess convergence of the cluster probabilities
  assess_convergence(model_fit, type = "cluster_probs")
  # Based on this, we set burnin = 2500
  # We now plot the posterior density of the scale parameters alpha in
  # each mixture:
  plot(model_fit, burnin = 2500, type = "alpha")
  # We can also compute the posterior density of the cluster probabilities
  plot(model_fit, burnin = 2500, type = "cluster_probs")
  # We can also plot the posterior cluster assignment. In this case, the assessors
  # are sorted according to their maximum a prior cluster estimate.
  plot(model_fit, burnin = 2500, type = "cluster_assignment")
  }





