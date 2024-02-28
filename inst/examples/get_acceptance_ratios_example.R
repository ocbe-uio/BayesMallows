library(parallel)
cl <- makeCluster(2)
set.seed(1)
mod <- compute_mallows(
  data = setup_rank_data(potato_visual),
  compute_options = set_compute_options(burnin = 200),
  cl = cl
)
stopCluster(cl)

# One ratio is returned for each cluster
get_acceptance_ratios(mod)
