# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
# model:
model_fit <- compute_mallows(potato_visual)

# Se the documentation to compute_mallows for how to assess the convergence of the algorithm
# Having chosen burin = 1000, we compute posterior intervals
model_fit$burnin <- 1000
# We then compute the CP consensus.
compute_cp_consensus(model_fit)

\dontrun{
  # CLUSTERWISE CP CONSENSUS
  # We can run a mixture of Mallows models, using the n_clusters argument
  # We use the sushi example data. See the documentation of compute_mallows for a more elaborate
  # example
  model_fit <- compute_mallows(sushi_rankings, n_clusters = 5)
  # Keeping the burnin at 1000, we can compute the CP consensus per cluster
  model_fit$burnin <- 1000
  cp_consensus_df <- compute_cp_consensus(model_fit)
  # Using dplyr::select and tidyr::cumprob we can now make a table
  # which shows the ranking in each cluster:
  library(dplyr)
  library(tidyr)
  cp_consensus_df %>%
    select(-cumprob) %>%
    spread(key = cluster, value = item)
}


