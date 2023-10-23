set.seed(994)

n_items <- dim(sushi_rankings)[2]
metric <- "footrule"

# Perform the resample-move SMC algorithm
smc_test <- smc_mallows_new_users(
  rankings = sushi_rankings[1:100, ], type = "complete",
  metric = metric, leap_size = floor(n_items / 5), N = 100, Time = 10,
  mcmc_kernel_app = 5, num_new_obs = 5,
  alpha_prop_sd = 0.5, lambda = 0.15, alpha_max = 1e3
)

# Plot rho
plot(smc_test, colnames = colnames(sushi_rankings), parameter = "rho")

# Plot alpha
plot(smc_test, parameter = "alpha")
