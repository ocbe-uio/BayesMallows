library(BayesMallows)
library(dplyr)
library(ggplot2)
model_fit1 <- compute_mallows(sushi_rankings, nmc = 1000)
model_fit2 <- compute_mallows(sushi_rankings, n_clusters = 2, nmc = 1000)

ggsave(filename = "./vignettes/vignette-data/alpha_convergence_1.png", plot = assess_convergence(model_fit1))
ggsave(filename = "./vignettes/vignette-data/alpha_convergence_2.png", plot = assess_convergence(model_fit2))

cp_consensus_sushi <- compute_cp_consensus(model_fit2, burnin = 200)
save(cp_consensus_sushi, file = "./vignettes/vignette-data/cp_consensus_sushi.RData")
