library(BayesMallows)
library(tidyverse)
set.seed(1)
leap_size <- 5
example_dataset <- potato_visual
metric <- "footrule"
n_users <- 12
n_items <- 20
test_dataset <- array(0, c(n_users, n_items, (n_items / 2 + 1)))
test_dataset[, , (n_items / 2 + 1)] <- potato_visual
tt <- 0
for (ii in (n_items - 1):(n_items / 2)) {
  tt <- tt + 1

  # set n_users line with one more NA
  example_dataset[example_dataset > ii] <- NA

  # set as new time stamp
  test_dataset[, , ((n_items / 2 + 1) - tt)] <- example_dataset
}

all.equal(as.numeric(test_dataset[,, 11]),
          as.numeric(potato_visual))

logz_list <- prepare_partition_function(metric = metric, n_items = n_items)

Time <- dim(test_dataset)[3]
N <- 5000
aug_method <- "pseudolikelihood"
metric <- "footrule"
smc_test_updated_partial <- smc_mallows_new_item_rank(
  n_items = n_items,
  R_obs = test_dataset,
  metric = metric,
  leap_size = leap_size, N = N,
  Time = Time,
  logz_estimate = NULL,
  cardinalities = logz_list$cardinalities,
  mcmc_kernel_app = 20,
  alpha_prop_sd = 0.5,
  lambda = 0.15,
  alpha_max = 1e6,
  aug_method = aug_method
)

bmm_ref <- compute_mallows(potato_visual, nmc = 20000)
bmm_ref$burnin <- 1000

smc_consensus <- as.numeric(str_extract(compute_consensus(smc_test_updated_partial)$item, "[:digit:]+"))
bmm_consensus <- as.numeric(str_extract(compute_consensus(bmm_ref)$item, "[:digit:]+"))

plot(smc_consensus, bmm_consensus)

# Number of switches to get from one to the other
rank_distance(matrix(smc_consensus, nrow = 1), bmm_consensus, metric = "ulam")

mean(smc_test_updated_partial$alpha_samples[, 11])
mean(bmm_ref$alpha$value[bmm_ref$alpha$iteration > 1000])

smc_ci <- compute_posterior_intervals(smc_test_updated_partial, parameter = "rho")
smc_ci[order(as.integer(str_extract(smc_ci$item, "[:digit:]+"))), ]

compute_posterior_intervals(bmm_ref, parameter = "rho")
