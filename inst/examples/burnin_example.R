set.seed(445)
mod <- compute_mallows(setup_rank_data(potato_visual))
assess_convergence(mod)
burnin(mod)
burnin(mod) <- 1500
burnin(mod)
plot(mod)
#'
models <- compute_mallows_mixtures(
  data = setup_rank_data(cluster_data),
  n_clusters = 1:3)
burnin(models)
burnin(models) <- 100
burnin(models)
burnin(models) <- c(100, 300, 200)
burnin(models)
