# Here is an example with the beach preference data.
# First, generate the transitive closure.
beach_tc <- generate_transitive_closure(beach_preferences)

# Next, generate an initial ranking.
beach_init_rank <- generate_initial_ranking(beach_tc)

# Then generate the constrain set used intervally by compute_mallows
constr <- generate_constraints(beach_tc, n_items = 15)

# Provide all these elements to compute_mallows
model_fit <- compute_mallows(rankings = beach_init_rank,
preferences = beach_tc, constraints = constr)

\dontrun{
  # The constraints can also be generated in parallel
  library(parallel)
  cl <- makeCluster(detectCores() - 1)
  constr <- generate_constraints(beach_tc, n_items = 15, cl = cl)
  stopCluster(cl)
}
