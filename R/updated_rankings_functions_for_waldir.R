library("BayesMallows")
library("testthat")

source("smc_mallows_new_users_new_item_rank_partial_uniform.R")
source("smc_mallows_new_users_new_item_rank_partial_pseudo.R")

# updated rankings code with mcmc and smcc input initialiasation, otherwise do the standard?

# @title SMC-Mallows new item rank updated
# @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user at each time step given an initial particle set obtained from MCMC or SMC. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.

# @param n_items Integer is the number of items in a ranking.
# @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time steps.
# @param metric A character string specifying the distance metric to use in the Bayesian Mallows Model. Available options are \code{"footrule"},
# \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and \code{"ulam"}.
# @param leap_size leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
# @param N Integer specifying the number of particles.
# @param Time Integer specifying the number of time steps in the SMC algorithm.
# @param logz_estimate Estimate of the partition function, computed with \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
# @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel.
# @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha.
# @param lambda Strictly positive numeric value specifying the rate parameter of the truncated exponential prior distribution of alpha.
# @param alpha_max  Maximum value of alpha in the truncated exponential prior distribution.
# @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random".
# @param alpha_samples_init A vector of size N by containing the initial particle set values of alpha.
# @param rho_samples_init 2D matrix of size N by n_items containing the initial particle set values of rho.
# @param aug_rankings_init 3D matrix of size n_assessors by n_items by N containing the initial particle set values of augmented rankings for R_obs.


# @return a 3d matrix containing: the samples of: rho, alpha and the augmented rankings, and the effective sample size at each iteration of the SMC algorithm.


smc_mallows_new_item_rank_updated = function(n_items, R_obs, metric, leap_size, N, Time, logz_estimate,
                                             mcmc_kernel_app, alpha_prop_sd, lambda, alpha_max, aug_method,
                                             alpha_samples_init, rho_samples_init, aug_rankings_init){


  ######################
  ## Initialise Phase
  ######################
  # Generate N initial samples of rho using the uniform prior
  rho_samples = array(0,  c(N, n_items, Time))
  rho_samples[,,1] = rho_samples_init

  alpha_samples = array(0, c(N, Time))
  alpha_samples[,1] = alpha_samples_init

  # store ESS
  ESS_vec = rep(0, times = Time)
  ESS_vec[1] = 1

  ######################
  ## Augment Rankings
  ######################
  num_ranks = dim(R_obs[,,1])[1]

  # each particle has its own set of augmented rankings
  aug_rankings = array(0, c(num_ranks, n_items, N))
  prev_aug_rankings = array(0, c(num_ranks, n_items, N))

  # augment incomplete ranks to initialise
  ranks = c(1:n_items)

  # augment rankings for proposal
  for(ii in 1:N){
    aug_rankings[,,ii] = aug_rankings_init[,,ii]
  }
  ## set the old augmentation as we will compare the t^th augmentation to the (t-1)^th augmentation
  prev_aug_rankings[,,] = aug_rankings[,,]


  #########################
  ## Loop for t=1,...,Time
  #########################
  # Here, we attempt the SMC sampler
  for (tt in 1:(Time-1)){

    #print( paste("iteration", tt, "out of" , Time))

    ###########################
    ## New Information
    ###########################
    # new observed item ranks from each user, need to update augmented rankings
    rho_samples[,,tt+1] = rho_samples[,,tt]
    alpha_samples[,tt+1] = alpha_samples[,tt]

    # total correction prob
    particle_correction_prob = rep(1,N)

    # iterate through each observed ranking and create new "corrected" augmented rankings
    for (ii in 1:N){
      # set t-1 generation to old as we sample for t new
      prev_aug_rankings[,,ii] = aug_rankings[,,ii]

      # make the correction
      for (jj in 1:num_ranks){
        if(aug_method == "random"){
          check_correction  = BayesMallows:::correction_kernel(current_ranking = aug_rankings[jj,,ii],
                                                               observed_ranking = R_obs[jj,,tt+1],
                                                               n_items = n_items)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          check_correction = BayesMallows:::correction_kernel_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                     observed_ranking = R_obs[jj,,tt+1],
                                                                     rho = rho_samples[ii,,tt+1] ,
                                                                     alpha =  alpha_samples[ii,tt+1],
                                                                     n_items = n_items,
                                                                     metric = metric)
        }else{

          print("Error: combined choice of metric and aug_method is incompatible")
          abort(message = "The value is TRUE, so the script must end here")

        }
        #print(check_correction)
        aug_rankings[jj,,ii] = check_correction$ranking

        particle_correction_prob[ii] = particle_correction_prob[ii] * check_correction$correction_prob
      }
    }

    ###########################
    # Re-weight
    ###########################
    #incremental weight for each particle, based on new observed rankings
    log_inc_wgt = rep(0, N)

    for (ii in 1:N){
      log_inc_wgt[ii] = BayesMallows::get_mallows_loglik(alpha = alpha_samples[ii,tt+1],
                                                         rho = rho_samples[ii,,tt+1],
                                                         n_items =n_items,
                                                         rankings = aug_rankings[,,ii],
                                                         metric = metric) -
                        BayesMallows::get_mallows_loglik(alpha = alpha_samples[ii,tt+1],
                                                         rho = rho_samples[ii,,tt+1],
                                                         n_items =n_items,
                                                         rankings = prev_aug_rankings[,,ii],
                                                         metric = metric) -
                        log(particle_correction_prob[ii])
    }

    # update weights
    maxw = max(log_inc_wgt)
    w = exp(log_inc_wgt-maxw)
    norm_wgt = w/sum(w)

    # store ESS
    ESS_vec[tt+1] = sum(norm_wgt)^2/sum(norm_wgt^2)

    ##############
    # Resample
    ##############
    index = sample(1:N, prob=norm_wgt, size=N, replace=T)
    rho_samples[,,tt+1] = rho_samples[index,,tt+1]
    alpha_samples[,tt+1] = alpha_samples[index,tt+1]
    aug_rankings[,,] = aug_rankings[,,index]

    ##################
    # Move step
    ##################
    for (ii in 1:N){

      for (kk in 1: mcmc_kernel_app) {
        rho_samples[ii,,tt+1] = BayesMallows::metropolis_hastings_rho(alpha = alpha_samples[ii,tt+1],
                                                                      n_items = n_items,
                                                                      rankings = aug_rankings[,,ii],
                                                                      metric = metric,
                                                                      rho = rho_samples[ii,,tt+1],
                                                                      leap_size = leap_size)

        # move once since alpha dist is easier to explore than rho dist
        alpha_samples[ii,tt+1] = BayesMallows::metropolis_hastings_alpha(alpha = alpha_samples[ii,tt+1],
                                                                         n_items = n_items,
                                                                         rankings = aug_rankings[,,ii],
                                                                         metric = metric,
                                                                         rho = rho_samples[ii,,tt+1],
                                                                         logz_estimate = logz_estimate,
                                                                         alpha_prop_sd = alpha_prop_sd,
                                                                         lambda = lambda,
                                                                         alpha_max = alpha_max)
      }

      for (jj in 1:num_ranks){
        if(aug_method == "random"){
          aug_rankings[jj,,ii] = BayesMallows::metropolis_hastings_aug_ranking(current_ranking = aug_rankings[jj,,ii],
                                                                               partial_ranking = R_obs[jj,,tt+1],
                                                                               alpha = alpha_samples[ii,tt+1],
                                                                               rho = rho_samples[ii,,tt+1],
                                                                               n_items = n_items,
                                                                               metric = metric)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = BayesMallows:::metropolis_hastings_aug_ranking_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                                       partial_ranking = R_obs[jj,,tt+1],
                                                                                       alpha = alpha_samples[ii,tt+1],
                                                                                       rho = rho_samples[ii,,tt+1],
                                                                                       n_items = n_items,
                                                                                       metric = metric)
        }

      }
    }
  }

  ############################
  # Post Processing
  ############################
  smc_list <- list("rho_samples" = rho_samples,
                   "alpha_samples" = alpha_samples,
                   "augmented_rankings" = aug_rankings,
                   "ESS" = ESS_vec )
  return(smc_list)
}



# @title SMC-Mallows new item rank updated alpha fixed
# @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user at each time step given an initial particle set obtained from MCMC or SMC. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.

# @param alpha numeric value of the scale parameter.
# @param n_items Integer is the number of items in a ranking.
# @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time steps.
# @param metric A character string specifying the distance metric to use in the Bayesian Mallows Model. Available options are \code{"footrule"},
# \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and \code{"ulam"}.
# @param leap_size leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
# @param N Integer specifying the number of particles.
# @param Time Integer specifying the number of time steps in the SMC algorithm.
# @param logz_estimate Estimate of the partition function, computed with \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
# @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel.
# @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random".
# @param rho_samples_init 2D matrix of size N by n_items containing the initial particle set values of rho.
# @param aug_rankings_init 3D matrix of size n_assessors by n_items by N containing the initial particle set values of augmented rankings for R_obs.


# @return a 3d matrix containing the samples of rho and the augmented rankings, and the effective sample size at each iteration of the SMC algorithm.

smc_mallows_new_item_rank_updated_alpha_fixed = function(alpha, n_items, R_obs, metric, leap_size, N, Time, logz_estimate,
                                                         mcmc_kernel_app, aug_method, rho_samples_init, aug_rankings_init){


  ######################
  ## Initialise Phase
  ######################
  # Generate N initial samples of rho using the uniform prior
  rho_samples = array(0,  c(N, n_items, Time))
  rho_samples[,,1] = rho_samples_init

  # store ESS
  ESS_vec = rep(0, times = Time)
  ESS_vec[1] = 1

  ######################
  ## Augment Rankings
  ######################
  num_ranks = dim(R_obs[,,1])[1]

  # each particle has its own set of augmented rankings
  aug_rankings = array(0, c(num_ranks, n_items, N))
  prev_aug_rankings = array(0, c(num_ranks, n_items, N))

  # augment incomplete ranks to initialise
  ranks = c(1:n_items)

  # augment rankings for proposal
  for(ii in 1:N){
    aug_rankings[,,ii] = aug_rankings_init[,,ii]
  }
  ## set the old augmentation as we will compare the t^th augmentation to the (t-1)^th augmentation
  prev_aug_rankings[,,] = aug_rankings[,,]


  #########################
  ## Loop for t=1,...,Time
  #########################
  # Here, we attempt the SMC sampler
  for (tt in 1:(Time-1)){

    #print( paste("iteration", tt, "out of" , Time))

    ###########################
    ## New Information
    ###########################
    # new observed item ranks from each user, need to update augmented rankings
    rho_samples[,,tt+1] = rho_samples[,,tt]
    #alpha_samples[,tt+1] = alpha_samples[,tt]

    # total correction prob
    particle_correction_prob = rep(1,N)

    # iterate through each observed ranking and create new "corrected" augmented rankings
    for (ii in 1:N){
      # set t-1 generation to old as we sample for t new
      prev_aug_rankings[,,ii] = aug_rankings[,,ii]

      # make the correction
      for (jj in 1:num_ranks){
        if(aug_method == "random"){
          check_correction  = BayesMallows:::correction_kernel(current_ranking = aug_rankings[jj,,ii],
                                                               observed_ranking = R_obs[jj,,tt+1],
                                                               n_items = n_items)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          check_correction = BayesMallows:::correction_kernel_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                     observed_ranking = R_obs[jj,,tt+1],
                                                                     rho = rho_samples[ii,,tt+1] ,
                                                                     alpha =  alpha,
                                                                     n_items = n_items,
                                                                     metric = metric)
        }else{

          print("Error: combined choice of metric and aug_method is incompatible")
          abort(message = "The value is TRUE, so the script must end here")

        }
        #print(check_correction)
        aug_rankings[jj,,ii] = check_correction$ranking

        particle_correction_prob[ii] = particle_correction_prob[ii] * check_correction$correction_prob
      }
    }

    ###########################
    # Re-weight
    ###########################
    #incremental weight for each particle, based on new observed rankings
    log_inc_wgt = rep(0, N)

    for (ii in 1:N){
      log_inc_wgt[ii] = BayesMallows::get_mallows_loglik(alpha = alpha, rho = rho_samples[ii,,tt+1],
                                                         n_items =n_items, rankings = aug_rankings[,,ii], metric = metric) -
        BayesMallows::get_mallows_loglik(alpha = alpha, rho = rho_samples[ii,,tt+1],
                                         n_items =n_items, rankings = prev_aug_rankings[,,ii], metric = metric) -
        log(particle_correction_prob[ii])
    }

    # update weights
    maxw = max(log_inc_wgt)
    w = exp(log_inc_wgt-maxw)
    norm_wgt = w/sum(w)

    # store ESS
    ESS_vec[tt+1] = sum(norm_wgt)^2/sum(norm_wgt^2)

    ##############
    # Resample
    ##############
    index = sample(1:N, prob=norm_wgt, size=N, replace=T)
    rho_samples[,,tt+1] = rho_samples[index,,tt+1]
    aug_rankings[,,] = aug_rankings[,,index]

    ##################
    # Move step
    ##################
    for (ii in 1:N){

      for (kk in 1: mcmc_kernel_app) {
        rho_samples[ii,,tt+1] = BayesMallows::metropolis_hastings_rho(alpha = alpha, n_items = n_items,
                                                                      rankings = aug_rankings[,,ii], metric = metric,
                                                                      rho = rho_samples[ii,,tt+1], leap_size = leap_size)

      }

      for (jj in 1:num_ranks){
        if(aug_method == "random"){
          aug_rankings[jj,,ii] = BayesMallows::metropolis_hastings_aug_ranking(current_ranking = aug_rankings[jj,,ii],
                                                                               partial_ranking = R_obs[jj,,tt+1],
                                                                               alpha = alpha,
                                                                               rho = rho_samples[ii,,tt+1],
                                                                               n_items = n_items,
                                                                               metric = metric)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = BayesMallows:::metropolis_hastings_aug_ranking_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                                       partial_ranking = R_obs[jj,,tt+1],
                                                                                       alpha = alpha,
                                                                                       rho = rho_samples[ii,,tt+1],
                                                                                       n_items = n_items,
                                                                                       metric = metric)
        }
      }
    }
  }

  ############################
  # Post Processing
  ############################
  smc_list <- list("rho_samples" = rho_samples,
                   "augmented_rankings" = aug_rankings,
                   "ESS" = ESS_vec)
  return(smc_list)
}


#####################
# Unit tests
#####################


# set-up

# new item rank for each user (fewer things)
example_dataset <- sushi_rankings[1:100,]
n_users <- dim(example_dataset)[1]
n_items <- dim(example_dataset)[2]
test_dataset <- array(0,  c(n_users, n_items, (n_items/2+ 1)) )
test_dataset[,,(n_items/2 + 1)] <- example_dataset
tt = 0
for (ii in (n_items-1):(n_items/2)){

  tt = tt + 1

  # set n_users line with one more NA
  example_dataset[example_dataset > ii] <- NA

  # set as new time stamp
  test_dataset[,,((n_items/2+1) - tt)] <- example_dataset
}


# Generate estimate of Z_n(alpha) ==============================================
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e4
degree <- 10

metric = "footrule"
leap_size = floor(n_items/5)

# Estimate the logarithm of the partition function of the Mallows rank model
# using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)



# check metric and aug_method error
expect_error(
  smc_mallows_new_item_rank_updated_alpha_fixed(alpha = 2, n_items = n_items,
                                                R_obs = test_dataset, metric = "cayley", leap_size = leap_size,
                                                N = N, Time = Time2, logz_estimate = logz_estimate,
                                                mcmc_kernel_app = mcmc_kernel_app, aug_method = "pseudolikelihood",
                                                rho_samples_init = smc_test_new_user_unif$rho_samples[,,Time+1],
                                                aug_rankings_init = smc_test_new_user_unif$aug_rankings)
)

expect_error(
  smc_mallows_new_item_rank_updated(n_items = n_items,
                                    R_obs = test_dataset, metric = "cayley", leap_size = leap_size,
                                    N = N, Time = Time2, logz_estimate = logz_estimate,
                                    mcmc_kernel_app = mcmc_kernel_app, alpha_prop_sd = 0.5,
                                    lambda = 0.1, alpha_max = 20, aug_method = "pseudolikelihood",
                                    alpha_samples_init = smc_test_new_user_unif$alpha_samples[,Time+1],
                                    rho_samples_init = smc_test_new_user_unif$rho_samples[,,Time+1],
                                    aug_rankings_init = smc_test_new_user_unif$aug_rankings)
)



# test with random sampler
N = 1000
mcmc_kernel_app = 5
num_new_obs = 10
Time = n_users/num_new_obs
sample_dataset = example_dataset

# run smc new user with uniform
set.seed(994)
smc_test_new_user_unif = smc_mallows_new_users_uniform(R_obs = sample_dataset, n_items = n_items, metric = metric,
                                                       leap_size = leap_size,
                                                       N = N, Time = Time, logz_estimate = logz_estimate,
                                                       mcmc_kernel_app = mcmc_kernel_app, num_new_obs = num_new_obs,
                                                       alpha_prop_sd = 0.5, lambda = 0.1,
                                                       alpha_max = 20)

# run smc updated rankings with alpha unknown
Time2 = dim(test_dataset)[3]
set.seed(994)
smc_test_updated_partial_unif1 = smc_mallows_new_item_rank_updated_alpha_fixed(alpha = 2, n_items = n_items,
                                                                   R_obs = test_dataset, metric = metric, leap_size = leap_size,
                                                                   N = N, Time = Time2, logz_estimate = logz_estimate,
                                                                   mcmc_kernel_app = mcmc_kernel_app, aug_method = "random",
                                                                   rho_samples_init = smc_test_new_user_unif$rho_samples[,,Time+1],
                                                                   aug_rankings_init = smc_test_new_user_unif$aug_rankings)
expect_is(smc_test_updated_partial_unif1, "list")
expect_length(smc_test_updated_partial_unif1, 3)
expect_equal(dim(smc_test_updated_partial_unif1$rho_samples), c(N, n_items, 6))
expect_length(smc_test_updated_partial_unif1$ESS, Time2)
expect_equal(dim(smc_test_updated_partial_unif1$augmented_rankings), c(n_users, n_items, N))

# run smc updated rankings with alpha unknown
Time2 = dim(test_dataset)[3]
set.seed(994)
smc_test_updated_partial_unif2 = smc_mallows_new_item_rank_updated(n_items = n_items,
                                                                   R_obs = test_dataset, metric = metric, leap_size = leap_size,
                                                                   N = N, Time = Time2, logz_estimate = logz_estimate,
                                                                   mcmc_kernel_app = mcmc_kernel_app, alpha_prop_sd = 0.5,
                                                                   lambda = 0.1, alpha_max = 20, aug_method = "random",
                                                                   alpha_samples_init = smc_test_new_user_unif$alpha_samples[,Time+1],
                                                                   rho_samples_init = smc_test_new_user_unif$rho_samples[,,Time+1],
                                                                   aug_rankings_init = smc_test_new_user_unif$aug_rankings)
expect_is(smc_test_updated_partial_unif2, "list")
expect_length(smc_test_updated_partial_unif2, 4)
expect_equal(dim(smc_test_updated_partial_unif2$rho_samples), c(N, n_items, 6))
expect_length(smc_test_updated_partial_unif2$ESS, Time2)
expect_equal(dim(smc_test_updated_partial_unif2$augmented_rankings), c(n_users, n_items, N))
expect_equal(dim(smc_test_updated_partial_unif2$alpha_samples), c(N, 6))


# test with pseudolikelihood


N = 1000
mcmc_kernel_app = 5
num_new_obs = 10
Time = n_users/num_new_obs
sample_dataset = example_dataset
set.seed(994)
smc_test_new_user_pseudo = smc_mallows_new_users_pseudo(R_obs = example_dataset, n_items = n_items, metric = metric,
                                                        leap_size = leap_size,
                                                        N = N, Time = Time, logz_estimate = logz_estimate,
                                                        mcmc_kernel_app = mcmc_kernel_app, num_new_obs = num_new_obs,
                                                        alpha_prop_sd = 0.5, lambda = 0.1,
                                                        alpha_max = 20)


set.seed(994)
smc_test_updated_partial_pseudo1 = smc_mallows_new_item_rank_updated_alpha_fixed(alpha = 2, n_items = n_items,
                                                                                 R_obs = test_dataset, metric = metric, leap_size = leap_size,
                                                                                 N = N, Time = Time2, logz_estimate = logz_estimate,
                                                                                 mcmc_kernel_app = mcmc_kernel_app,  aug_method = "pseudolikelihood",
                                                                                 rho_samples_init = smc_test_new_user_pseudo$rho_samples[,,Time+1],
                                                                                 aug_rankings_init = smc_test_new_user_pseudo$aug_rankings)

expect_is(smc_test_updated_partial_pseudo1, "list")
expect_length(smc_test_updated_partial_pseudo1, 3)
expect_equal(dim(smc_test_updated_partial_pseudo1$rho_samples), c(N, n_items, 6))
expect_length(smc_test_updated_partial_pseudo1$ESS, Time2)
expect_equal(dim(smc_test_updated_partial_pseudo1$augmented_rankings), c(n_users, n_items, N))


set.seed(994)
smc_test_updated_partial_pseudo2 = smc_mallows_new_item_rank_updated(n_items = n_items,
                                                                   R_obs = test_dataset, metric = metric, leap_size = leap_size,
                                                                   N = N, Time = Time2, logz_estimate = logz_estimate,
                                                                   mcmc_kernel_app = mcmc_kernel_app, alpha_prop_sd = 0.5,
                                                                   lambda = 0.1, alpha_max = 20, aug_method = "pseudolikelihood",
                                                                   alpha_samples_init = smc_test_new_user_unif$alpha_samples[,Time+1],
                                                                   rho_samples_init = smc_test_new_user_unif$rho_samples[,,Time+1],
                                                                   aug_rankings_init = smc_test_new_user_unif$aug_rankings)

expect_is(smc_test_updated_partial_pseudo2, "list")
expect_length(smc_test_updated_partial_pseudo2, 4)
expect_equal(dim(smc_test_updated_partial_pseudo2$rho_samples), c(N, n_items, 6))
expect_length(smc_test_updated_partial_pseudo2$ESS, Time2)
expect_equal(dim(smc_test_updated_partial_pseudo2$augmented_rankings), c(n_users, n_items, N))
expect_equal(dim(smc_test_updated_partial_pseudo2$alpha_samples), c(N, 6))
