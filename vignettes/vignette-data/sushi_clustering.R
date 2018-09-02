library(BayesMallows)
library(dplyr)
library(ggplot2)

model_fit1 <- compute_mallows(sushi_rankings, nmc = 1000, include_wcd = TRUE)
model_fit2 <- compute_mallows(sushi_rankings, n_clusters = 2, nmc = 1000)

ggsave(filename = "./vignettes/vignette-data/alpha_convergence_2.png",
       plot = assess_convergence(model_fit2), dpi = 100)

ggsave(filename = "./vignettes/vignette-data/cluster_probs_convergence2.png",
       plot = assess_convergence(model_fit2, type = "cluster_probs"), dpi = 100)

cp_consensus_sushi <- compute_cp_consensus(model_fit2, burnin = 200)

n_clusters <- seq(from = 1, to = 10)
models <- compute_mallows_mixtures(n_clusters = n_clusters,
                                   rankings = sushi_rankings,
                                   nmc = 6000)

ggsave(filename = "./vignettes/vignette-data/sushi_elbow_10.png",
       plot = plot_elbow(models, burnin = 1000), dpi = 100)

cluster_assignment <- assign_cluster(models[[5]], burnin = 400, soft = FALSE)
cluster_assignment <- head(cluster_assignment)

save(cp_consensus_sushi, cluster_assignment,
     file = "./vignettes/vignette-data/sushi_clustering_data.RData")

ggsave(filename = "./vignettes/vignette-data/cluster_assignment_plot_10.png",
       plot = plot(models[[5]], burnin = 1000, type = "cluster_assignment"),
       dpi = 100)

