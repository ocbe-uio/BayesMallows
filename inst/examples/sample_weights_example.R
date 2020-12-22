# The first example uses full rankings in the potato_visual dataset, but we assume
# that each row in the data corresponds to between 100 and 500 assessors.
set.seed(1234)
# We start by generating random sample weights
weights <- sample(x = seq(from = 100L, to = 500L, by = 1L),
                  size = nrow(potato_visual), replace = TRUE)
# We also create a set of repeated indices, used to extend the matrix rows
repeated_indices <- unlist(purrr::map2(1:nrow(potato_visual), weights, ~ rep(.x, each = .y)))
# The potato_repeated matrix consists of all rows repeated corresponding to
# the number of assessors in the weights vector. This is how a large dataset
# would look like without using the weights argument
potato_repeated <- potato_visual[repeated_indices, ]

# We now first compute the Mallows model using weights
# This takes about 0.2 seconds
system.time({
  m_weights <- compute_mallows(rankings = potato_visual, weights = weights, nmc = 10000)
})
# Next we use the full ranking matrix
# This takes about 11.3 seconds, about 50 times longer!
\dontrun{
system.time({
  m_rep <- compute_mallows(rankings = potato_repeated, nmc = 10000)
})

  # We set the burnin to 2000 for both
  m_weights$burnin <- 2000
  m_rep$burnin <- 2000

  # Note that the MCMC algorithms did not run with the same
  # random number seeds in these two experiments, but still
  # the posterior distributions look similar
  plot(m_weights, burnin = 2000, "alpha")
  plot(m_rep, burnin = 2000, "alpha")

  plot(m_weights, burnin = 2000, "rho", items = 1:4)
  plot(m_rep, burnin = 2000, "rho", items = 1:4)
}

# Next we repeated the exercise with the pairwise preference data
# in the beach dataset
