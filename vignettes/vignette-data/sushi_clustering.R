library(BayesMallows)
library(dplyr)
library(ggplot2)
model_fit1 <- compute_mallows(sushi_rankings, nmc = 1000)
model_fit2 <- compute_mallows(sushi_rankings, n_clusters = 2, nmc = 1000)

ggsave(filename = "./vignettes/vignette-data/alpha_convergence_1.png", plot = assess_convergence(model_fit1))
ggsave(filename = "./vignettes/vignette-data/alpha_convergence_2.png", plot = assess_convergence(model_fit2))

cp_consensus_sushi <- compute_cp_consensus(model_fit2, burnin = 200)

clus_dis1 <- compute_within_cluster_distance(model_fit1, rankings = sushi_rankings, burnin = 200)
clus_dis2 <- compute_within_cluster_distance(model_fit2, rankings = sushi_rankings, burnin = 200)

ggsave(filename = "./vignettes/vignette-data/sushi_elbow.png",
       plot_elbow(model_fit1, model_fit2, rankings = sushi_rankings, burnin = 200))

save(cp_consensus_sushi, clus_dis1, clus_dis2, file = "./vignettes/vignette-data/sushi_clustering_data.RData")
