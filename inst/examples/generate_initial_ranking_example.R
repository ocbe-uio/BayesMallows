# The example dataset beach_preferences contains pairwise prefences of beach.
# We must first generate the transitive closure
beach_tc <- generate_transitive_closure(beach_preferences)

# Next, we generate an initial ranking
beach_init <- generate_initial_ranking(beach_tc)

# Look at the first few rows:
head(beach_init)

# We can add more informative column names in order
# to get nicer posterior plots:
colnames(beach_init) <- paste("Beach", seq(from = 1, to = ncol(beach_init), by = 1))
head(beach_init)

\dontrun{
  # We now give beach_init and beach_tc to compute_mallows:
  model_fit <- compute_mallows(rankings = beach_init, preferences = beach_tc, nmc = 2000)

  # We can study the acceptance rate of the augmented rankings
  assess_convergence(model_fit, type = "augmentation")

  # We can also study the posterior distribution of the consensus rank of each beach
  plot(model_fit, type = "rho", burnin <- 500, items = 1:15)
}
