# test script

# smc new user and new item rank combined
require(BayesMallows)
require(tibble)
require(Rcpp)
require(ggplot2)
require(plotrix)
require(purrr)
require(crayon)
require(utf8)
require(fields)
require(tidyr)
require(dplyr)
require(gridExtra)
require(prodlim)

source("leap_and_shift_probs.R")
source("get_mallows_loglik.R")
source("metropolis_hastings_rho.R")
source("metropolis_hastings_alpha.R")
# source("post_processing_functions.R")
source("uniform_functions.R")
source("pseudolikelihood_functions.R")
source("smc_mallows_new_item_rank.R")
# source("calculate_l1_distance_true_posterior.R")



# a simpler example to test
set.seed(101)
load("data_10_6.Rda")

Time <- dim(sample_dataset)[3]

# General
n_items <- dim(sample_dataset)[2] # Number of items
rho_0 <- seq(from = 1, to = n_items, by = 1) # 'true' consensus ranking
alpha_0 <- 2 # fixed/ 'true' scale parameter
leap_size <- floor(n_items / 5)
metric <- "footrule"

# Generate estimate of Z_n(alpha)
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e4
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)


mcmc_kernel_app <- 5
N <- 1000
alpha_prop_sd <- 0.5
lambda <- 0.15
alpha_max <- 1e6

source("smc_mallows_new_item_rank.R")

# check it will produce the wrong metric and aug_method error
smc_unif_alpha_fixed_unif <- smc_mallows_new_item_rank_alpha_fixed(
  alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
  metric = "cayley", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "pseudolikelihood"
)

smc_unif <- smc_mallows_new_item_rank(
  n_items = n_items, R_obs = sample_dataset,
  metric = "cayley", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "pseudolikelihood"
)

# check it runs with unif kernel
smc_unif_alpha_fixed_unif <- smc_mallows_new_item_rank_alpha_fixed(
  alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
  metric = "footrule", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "random"
)

smc_unif <- smc_mallows_new_item_rank(
  n_items = n_items, R_obs = sample_dataset,
  metric = "footrule", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "random"
)

# check it runs with pseudo kernel
smc_unif_alpha_fixed_unif <- smc_mallows_new_item_rank_alpha_fixed(
  alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
  metric = "footrule", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "pseudolikelihood"
)

smc_unif <- smc_mallows_new_item_rank(
  n_items = n_items, R_obs = sample_dataset,
  metric = "footrule", leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
  alpha_prop_sd = alpha_prop_sd, lambda = lambda,
  alpha_max = alpha_max, aug_method = "pseudolikelihood"
)




###################################################################
# comparing and testing code - a much more complex example to test
###################################################################
set.seed(101)
items <- c(6, 8)
n_rankings <- c(10, 20, 30)
particles <- c(1000, 5000, 10000)

for (ii in 1:2) {
  n_items <- items[ii]
  print(n_items)

  if (n_items == 6) {
    rho0 <- c(1, 4, 2, 5, 3, 6)
  } else {
    rho0 <- c(1, 5, 2, 6, 3, 7, 4, 8)
  }

  for (jj in 1:3) {
    rankings <- n_rankings[jj]
    print(rankings)
    load(paste0("data_", rankings, "_", n_items, ".Rda"))

    n_samples <- dim(sample_dataset)[1]
    n_items <- dim(sample_dataset)[2]
    Time <- dim(sample_dataset)[3]

    # General
    rho_0 <- seq(from = 1, to = n_items, by = 1) # true consensus ranking
    alpha_0 <- 2 # true scale parameter
    alpha <- alpha_0
    leap_size <- min(1, floor(n_items / 5)) # Leap size
    metric <- "footrule"

    # Estimate the logarithm of the partition function of the Mallows rank model.
    # We create a grid of alpha values from 0 to 10
    alpha_vector <- seq(from = 0, to = 20, by = 0.1)
    iter <- 1e4
    degree <- 10


    # Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
    logz_estimate <- estimate_partition_function(
      method = "importance_sampling",
      alpha_vector = alpha_vector,
      n_items = n_items, metric = metric,
      nmc = iter, degree = degree
    )

    mcmc_kernel_app <- 5

    alpha_prop_sd <- 0.5
    lambda <- 0.15
    alpha_max <- 1e6

    for (kk in 1:3) {
      N <- particles[kk]
      print(N)

      set.seed(101)


      #######################################
      # uniform new item rank alpha is fixed
      #######################################
      print("unif new item rank alpha fixed")
      smc_unif_alpha_fixed <- smc_mallows_new_item_rank_alpha_fixed(
        alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
        metric = metric, leap_size = leap_size, N = N, Time = Time,
        logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
        alpha_prop_sd = alpha_prop_sd, lambda = lambda,
        alpha_max = alpha_max, aug_method = "random"
      )
      save(smc_unif_alpha_fixed, file = paste0("output_", n_samples, "_", n_items, "_", N, "unif_alpha_fixed.Rda"))



      ###############################
      ## alpha unknown
      ##############################
      print("uniform new item rank alpha unknown")
      smc_unif_both <- smc_mallows_new_item_rank(
        n_items = n_items, R_obs = sample_dataset,
        metric = metric, leap_size = leap_size, N = N, Time = Time,
        logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
        alpha_prop_sd = alpha_prop_sd, lambda = lambda,
        alpha_max = alpha_max, aug_method = "random"
      )
      save(smc_unif_both, file = paste0("output_", n_samples, "_", n_items, "_", N, "unif_both.Rda"))


      ###################################################
      # pseudo new item rank alpha fixed
      ####################################################
      print("pseudo new item rank but alpha fixed")
      smc_pseudo_alpha_fixed <- smc_mallows_new_item_rank_alpha_fixed(
        alpha = alpha_0, n_items = n_items, R_obs = sample_dataset,
        metric = metric, leap_size = leap_size, N = N, Time = Time,
        logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
        alpha_prop_sd = alpha_prop_sd, lambda = lambda,
        alpha_max = alpha_max, aug_method = "pseudolikelihood"
      )
      save(smc_pseudo_alpha_fixed, file = paste0("output_", n_samples, "_", n_items, "_", N, "pseudo_alpha_fixed.Rda"))


      #######################################################################
      # pseudo alpha unknown
      #######################################################################
      print("pseudo alpha unknown")
      smc_pseudo_both <- smc_mallows_new_item_rank(
        n_items = n_items, R_obs = sample_dataset,
        metric = metric, leap_size = leap_size, N = N, Time = Time,
        logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_kernel_app,
        alpha_prop_sd = alpha_prop_sd, lambda = lambda,
        alpha_max = alpha_max, aug_method = "pseudolikelihood"
      )
      save(smc_pseudo_both, file = paste0("output_", n_samples, "_", n_items, "_", N, "pseudo_both.Rda"))
    }
    rm(sample_dataset)
  }
}
