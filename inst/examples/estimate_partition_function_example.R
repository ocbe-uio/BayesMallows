# Let us estimate logZ(alpha) for 5 items with Spearman distance
# We create a grid of alpha values from 0 to 10
alpha_vector <- seq(from = 0.1, to = 10, by = 0.1)
# We start with 1e3 Monte Carlo samples
\dontrun{
# We take 1e4 Monte Carlo samples, and fit log(Z) with a 10-degree polynomial
estimate <- estimate_partition_function(method = "importance_sampling",
                                         alpha_vector = alpha_vector,
                                         n_items = 5,
                                         metric = "spearman",
                                         nmc = 1e3,
                                         degree = 10)
}
