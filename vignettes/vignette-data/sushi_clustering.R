library(BayesMallows)
library(dplyr)
library(ggplot2)

model_fit1 <- compute_mallows(sushi_rankings, nmc = 1000, include_wcd = TRUE)
model_fit2 <- compute_mallows(sushi_rankings, n_clusters = 2, nmc = 1000)

ggsave(filename = "./vignettes/vignette-data/alpha_convergence_1.png",
       plot = assess_convergence(model_fit1))

ggsave(filename = "./vignettes/vignette-data/alpha_convergence_2.png",
       plot = assess_convergence(model_fit2))

cp_consensus_sushi <- compute_cp_consensus(model_fit2, burnin = 200)

ggsave(filename = "./vignettes/vignette-data/sushi_elbow.png",
       plot_elbow(model_fit1, model_fit2, rankings = sushi_rankings, burnin = 200))

library(purrr)
n_clusters <- seq(from = 1, to = 10)
models <- map(n_clusters, ~ compute_mallows(rankings = sushi_rankings, nmc = 1000,
                                            n_clusters = .x, include_wcd = TRUE))

ggsave(filename = "./vignettes/vignette-data/sushi_elbow_10.png",
       plot = plot_elbow(models, burnin = 200))

save(cp_consensus_sushi, file = "./vignettes/vignette-data/sushi_clustering_data.RData")
