# The example dataset beach_preferences contains pairwise preferences of beach.
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

# It is also possible to pick a random sample among all topological sorts.
# This requires first enumerating all possible sorts, and might hence be computationally
# demanding. Here is an example, where we reduce the data considerable to speed up computation.
small_tc <- beach_tc[beach_tc$assessor %in% 1:6 & beach_tc$bottom_item %in% 1:4 & beach_tc$top_item %in% 1:4, ]
set.seed(123)
init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
# Look at the initial rankings generated
init_small

\dontrun{
  # We now give beach_init and beach_tc to compute_mallows. We tell compute_mallows
  # to save the augmented data, in order to study the convergence.
  model_fit <- compute_mallows(rankings = beach_init,
                               preferences = beach_tc,
                               nmc = 2000,
                               save_aug = TRUE)

  # We can study the acceptance rate of the augmented rankings
  assess_convergence(model_fit, parameter = "Rtilde")

  # We can also study the posterior distribution of the consensus rank of each beach
  model_fit$burnin <- 500
  plot(model_fit, parameter = "rho", items = 1:15)
}

\dontrun{
  # The computations can also be done in parallel
  library(parallel)
  cl <- makeCluster(detectCores() - 1)
  beach_tc <- generate_transitive_closure(beach_preferences, cl = cl)
  beach_init <- generate_initial_ranking(beach_tc, cl = cl)
  stopCluster(cl)
}
