# Generate basic elements
data <- sushi_rankings[1:100, ]
n_items <- ncol(sushi_rankings)
metric <- "footrule"
num_new_obs <- 10
logz_estimate <- estimate_partition_function(
	method = "importance_sampling",
	alpha_vector = seq(from = 0, to = 15, by = 0.1),
	n_items = n_items, metric = metric, nmc = 1e2, degree = 10
)

# Calculating rho and alpha samples
samples <- smc_mallows_new_users_complete(
	R_obs = data, n_items = n_items, metric = metric,
	leap_size = floor(n_items / 5), N = 100, Time = nrow(data) / num_new_obs,
	mcmc_kernel_app = 5, logz_estimate = logz_estimate,
	alpha_prop_sd = 0.1, lambda = 0.001, alpha_max = 1e6,
	num_new_obs = num_new_obs, verbose = TRUE
)

# Studying the structure of the output
str(samples)
