rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6
rankings <- sample_mallows(
  rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
)
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e2
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model
# using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling", alpha_vector = alpha_vector,
  n_items = n_items, metric = "footrule", nmc = iter, degree = degree
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate, alpha_prop_sd = 0.5,
  lambda = 0.1, alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate,
  alpha_prop_sd = 0.15, lambda = 0.1, alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate,
  alpha_prop_sd = 0.5, lambda = 0.15, alpha_max = 20, metric
)

metropolis_hastings_alpha(
  alpha, n_items, rankings, rho, logz_estimate,
  alpha_prop_sd = 0.15, lambda = 0.15, alpha_max = 20, metric
)
