model_fit <- compute_mallows(setup_rank_data(potato_visual))
burnin(model_fit) <- 1000

# By default, the scale parameter "alpha" is plotted
plot(model_fit)
# We can also plot the latent rankings "rho"
plot(model_fit, parameter = "rho")
# By default, a random subset of 5 items are plotted
# Specify which items to plot in the items argument.
plot(model_fit, parameter = "rho",
     items = c(2, 4, 6, 9, 10, 20))
# When the ranking matrix has column names, we can also
# specify these in the items argument.
# In this case, we have the following names:
colnames(potato_visual)
# We can therefore get the same plot with the following call:
plot(model_fit, parameter = "rho",
     items = c("P2", "P4", "P6", "P9", "P10", "P20"))

\dontrun{
  # Plots of mixture parameters:
  model_fit <- compute_mallows(
    setup_rank_data(sushi_rankings),
    model_options = set_model_options(n_clusters = 5))
  burnin(model_fit) <- 1000
  # Posterior distributions of the cluster probabilities
  plot(model_fit, parameter = "cluster_probs")
  # Cluster assignment plot. Color shows the probability of belonging to each
  # cluster.
  plot(model_fit, parameter = "cluster_assignment")
}




