# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
# model:
model_fit <- compute_mallows(potato_visual)

# Se the documentation to compute_mallows for how to assess the convergence
# of the algorithm
# We set the burnin = 1000
model_fit$burnin <- 1000
# By default, the scale parameter "alpha" is plotted
plot(model_fit)
\dontrun{
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
  }

\dontrun{
  # Plots of mixture parameters:
  # We can run a mixture of Mallows models, using the n_clusters argument
  # We use the sushi example data. See the documentation of compute_mallows for a more elaborate
  # example
  model_fit <- compute_mallows(sushi_rankings, n_clusters = 5, save_clus = TRUE)
  model_fit$burnin <- 1000
  # We can then plot the posterior distributions of the cluster probabilities
  plot(model_fit, parameter = "cluster_probs")
  # We can also get a cluster assignment plot, showing the assessors along the horizontal
  # axis and the clusters along the vertical axis. The color show the probability
  # of belonging to each clusters. The assessors are sorted along the horizontal
  # axis according to their maximum a posterior cluster assignment. This plot
  # illustrates the posterior uncertainty in cluster assignments.
  plot(model_fit, parameter = "cluster_assignment")
  # See also ?assign_cluster for a function which returns the cluster assignment
  # back in a dataframe.
}




