rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6
rankings <- sample_mallows(
  rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
)

# Copmute exact partition function
cardinalities <- prepare_partition_function(metric = metric,
                                            n_items = n_items)$cardinalities

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate = NULL,
  cardinalities = cardinalities, alpha_prop_sd = 0.5,
  lambda = 0.1, alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate = NULL,
  cardinalities = cardinalities, alpha_prop_sd = 0.15, lambda = 0.1,
  alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate = NULL,
  cardinalities = cardinalities, alpha_prop_sd = 0.5, lambda = 0.15,
  alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate = NULL,
  cardinalities = cardinalities, alpha_prop_sd = 0.15, lambda = 0.15,
  alpha_max = 20, metric
)
