# This file is here to precompute stuff in the vignette that take too long
library(BayesMallows)
library(ggplot2)

## Clustering of sushi data
sushi_clustering1 <- compute_mallows(sushi_rankings)
ggsave(filename = "./vignettes/vignette_data/sushi_convergence_1.png",
       plot = assess_convergence(sushi_clustering1))


sushi_clustering2 <- compute_mallows(sushi_rankings, n_clusters = 2)
ggsave(filename = "./vignettes/vignette_data/sushi_convergence_2.png",
       plot = assess_convergence(sushi_clustering2))
