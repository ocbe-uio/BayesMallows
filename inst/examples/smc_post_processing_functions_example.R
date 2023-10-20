# Basic elements
data <- sushi_rankings[1:100, ]
n_items <- ncol(data)
leap_size <- floor(n_items / 5)
metric <- "footrule"
Time <- 20
N <- 100

# Performing SMC
smc_test <- smc_mallows_new_users(
  R_obs = data, type = "complete", n_items = n_items,
  metric = metric, leap_size = leap_size,
  N = N, Time = Time,
  mcmc_kernel_app = 5,
  num_new_obs = 5,
  alpha_prop_sd = 0.5,
  lambda = 0.15,
  alpha_max = 1e6
)

compute_posterior_intervals(smc_test, parameter = "rho")

compute_consensus(model_fit = smc_test, type = "CP")
compute_consensus(model_fit = smc_test, type = "MAP")

compute_posterior_intervals(smc_test, parameter = "alpha")
