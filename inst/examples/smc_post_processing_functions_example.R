# Basic elements
data <- sushi_rankings[1:100, ]
n_items <- ncol(data)
leap_size <- floor(n_items / 5)
metric <- "footrule"
alpha_vector <- seq(from = 0, to = 15, by = 0.1)
iter <- 100
degree <- 10
N <- 100
Time <- 20

# Estimating the logarithm of the partition function of the Mallows rank model
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)

# Performing SMC
smc_test <- smc_mallows_new_users(
  R_obs = data, type = "complete", n_items = n_items,
  metric = metric, leap_size = leap_size,
  N = N, Time = Time,
  logz_estimate = logz_estimate,
  mcmc_kernel_app = 5,
  num_new_obs = 5,
  alpha_prop_sd = 0.5,
  lambda = 0.15,
  alpha_max = 1e6
)

test_sample_rho <- smc_test$rho_samples[, , Time + 1]
compute_posterior_intervals_rho(
  output = test_sample_rho,
  nmc = N, burnin = 0,
  verbose = FALSE
)

compute_rho_consensus(
  output = test_sample_rho, nmc = N,
  burnin = 0, C = 1, type = "CP",
  verbose = FALSE
)

test_sample_alpha <- smc_test$alpha_samples[, Time + 1]
compute_posterior_intervals_alpha(
  output = test_sample_alpha,
  nmc = N, burnin = 0, verbose = FALSE
)
