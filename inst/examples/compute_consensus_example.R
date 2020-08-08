# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
# model:
model_fit <- compute_mallows(potato_visual)

# Se the documentation to compute_mallows for how to assess the convergence of the algorithm
# Having chosen burin = 1000, we compute posterior intervals
model_fit$burnin <- 1000
# We then compute the CP consensus.
compute_consensus(model_fit, type = "CP")
# And we compute the MAP consensus
compute_consensus(model_fit, type = "MAP")

\dontrun{
  # CLUSTERWISE CONSENSUS
  # We can run a mixture of Mallows models, using the n_clusters argument
  # We use the sushi example data. See the documentation of compute_mallows for a more elaborate
  # example
  model_fit <- compute_mallows(sushi_rankings, n_clusters = 5)
  # Keeping the burnin at 1000, we can compute the consensus ranking per cluster
  model_fit$burnin <- 1000
  cp_consensus_df <- compute_consensus(model_fit, type = "CP")
  # Using dplyr::select and tidyr::cumprob we can now make a table
  # which shows the ranking in each cluster:
  library(dplyr)
  library(tidyr)
  cp_consensus_df %>%
    select(-cumprob) %>%
    spread(key = cluster, value = item)
}

\dontrun{
  # MAP CONSENSUS FOR PAIRWISE PREFENCE DATA
  # We use the example dataset with beach preferences.
  model_fit <- compute_mallows(preferences = beach_preferences)
  # We set burnin = 1000
  model_fit$burnin <- 1000
  # We now compute the MAP consensus
  map_consensus_df <- compute_consensus(model_fit, type = "MAP")
}

\dontrun{
  # CP CONSENSUS FOR AUGMENTED RANKINGS
  # We use the example dataset with beach preferences.
  model_fit <- compute_mallows(preferences = beach_preferences, save_aug = TRUE,
                               aug_thinning = 2, seed = 123L)
  # We set burnin = 1000
  model_fit$burnin <- 1000
  # We now compute the CP consensus of augmented ranks for assessors 1 and 3
  cp_consensus_df <- compute_consensus(model_fit, type = "CP",
                                       parameter = "Rtilde", assessors = c(1L, 3L))
  # We can also compute the MAP consensus for assessor 2
  map_consensus_df <- compute_consensus(model_fit, type = "MAP",
                                        parameter = "Rtilde", assessors = 2L)

  # Caution!
  # With very sparse data or with too few iterations, there may be ties in the MAP consensus
  # This is illustrated below for the case of only 5 post-burnin iterations. Two MAP rankings are
  # equally likely in this case (and for this seed).
  model_fit <- compute_mallows(preferences = beach_preferences, nmc = 1005,
                               save_aug = TRUE, aug_thinning = 1, seed = 123L)
  model_fit$burnin <- 1000
  compute_consensus(model_fit, type = "MAP", parameter = "Rtilde", assessors = 2L)
}
