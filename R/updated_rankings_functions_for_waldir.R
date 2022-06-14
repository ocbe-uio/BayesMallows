# updated rankings code with mcmc and smcc input initialiasation, otherwise do the standard?

#' @title SMC-Mallows new item rank updated
#' @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user at each time step given an initial particle set obtained from MCMC or SMC. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.
#' @param n_items Integer is the number of items in a ranking.
#' @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time steps.
#' @param metric A character string specifying the distance metric to use in the Bayesian Mallows Model. Available options are \code{"footrule"},
# \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and \code{"ulam"}.
#' @param leap_size leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
#' @param N Integer specifying the number of particles.
#' @param Time Integer specifying the number of time steps in the SMC algorithm.
#' @param logz_estimate Estimate of the partition function, computed with \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
#' @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel.
#' @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha.
#' @param lambda Strictly positive numeric value specifying the rate parameter of the truncated exponential prior distribution of alpha.
#' @param alpha_max  Maximum value of alpha in the truncated exponential prior distribution.
#' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random".
#' @param alpha_samples_init A vector of size N by containing the initial particle set values of alpha.
#' @param rho_samples_init 2D matrix of size N by n_items containing the initial particle set values of rho.
#' @param aug_rankings_init 3D matrix of size n_assessors by n_items by N containing the initial particle set values of augmented rankings for R_obs.
#' @return a 3d matrix containing: the samples of: rho, alpha and the augmented rankings, and the effective sample size at each iteration of the SMC algorithm.
#' @export
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
          check_correction  = correction_kernel(current_ranking = aug_rankings[jj,,ii],
                                                               observed_ranking = R_obs[jj,,tt+1],
                                                               n_items = n_items)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          check_correction = correction_kernel_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                     observed_ranking = R_obs[jj,,tt+1],
                                                                     rho = rho_samples[ii,,tt+1] ,
                                                                     alpha =  alpha_samples[ii,tt+1],
                                                                     n_items = n_items,
                                                                     metric = metric)
        }else{
          stop("Combined choice of metric and aug_method is incompatible. The value is TRUE, so the script must end here.")
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
      log_inc_wgt[ii] = get_exponent_sum(alpha = alpha_samples[ii,tt+1],
                                                         rho = rho_samples[ii,,tt+1],
                                                         n_items =n_items,
                                                         rankings = aug_rankings[,,ii],
                                                         metric = metric) -
                        get_exponent_sum(alpha = alpha_samples[ii,tt+1],
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
        rho_samples[ii,,tt+1] = metropolis_hastings_rho(alpha = alpha_samples[ii,tt+1],
                                                                      n_items = n_items,
                                                                      rankings = aug_rankings[,,ii],
                                                                      metric = metric,
                                                                      rho = rho_samples[ii,,tt+1],
                                                                      leap_size = leap_size)

        # move once since alpha dist is easier to explore than rho dist
        alpha_samples[ii,tt+1] = metropolis_hastings_alpha(alpha = alpha_samples[ii,tt+1],
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
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking(current_ranking = aug_rankings[jj,,ii],
                                                                               partial_ranking = R_obs[jj,,tt+1],
                                                                               alpha = alpha_samples[ii,tt+1],
                                                                               rho = rho_samples[ii,,tt+1],
                                                                               n_items = n_items,
                                                                               metric = metric)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking_pseudo(current_ranking = aug_rankings[jj,,ii],
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



#' @title SMC-Mallows new item rank updated alpha fixed
#' @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user at each time step given an initial particle set obtained from MCMC or SMC. Each correction and augmentation is done by filling in the missing item ranks using pseudolikelihood augmentation.

#' @param alpha numeric value of the scale parameter.
#' @param n_items Integer is the number of items in a ranking.
#' @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time steps.
#' @param metric A character string specifying the distance metric to use in the Bayesian Mallows Model. Available options are \code{"footrule"},
# \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and \code{"ulam"}.
#' @param leap_size leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
#' @param N Integer specifying the number of particles.
#' @param Time Integer specifying the number of time steps in the SMC algorithm.
#' @param logz_estimate Estimate of the partition function, computed with \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
#' @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel.
#' @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random".
#' @param rho_samples_init 2D matrix of size N by n_items containing the initial particle set values of rho.
#' @param aug_rankings_init 3D matrix of size n_assessors by n_items by N containing the initial particle set values of augmented rankings for R_obs.
#' @return a 3d matrix containing the samples of rho and the augmented rankings, and the effective sample size at each iteration of the SMC algorithm.
#' @export
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
          check_correction  = correction_kernel(current_ranking = aug_rankings[jj,,ii],
                                                               observed_ranking = R_obs[jj,,tt+1],
                                                               n_items = n_items)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          check_correction = correction_kernel_pseudo(current_ranking = aug_rankings[jj,,ii],
                                                                     observed_ranking = R_obs[jj,,tt+1],
                                                                     rho = rho_samples[ii,,tt+1] ,
                                                                     alpha =  alpha,
                                                                     n_items = n_items,
                                                                     metric = metric)
        }else{

          stop("Combined choice of metric and aug_method is incompatible. The value is TRUE, so the script must end here.")

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
      log_inc_wgt[ii] = get_exponent_sum(alpha = alpha, rho = rho_samples[ii,,tt+1],
                                                         n_items =n_items, rankings = aug_rankings[,,ii], metric = metric) -
        get_exponent_sum(alpha = alpha, rho = rho_samples[ii,,tt+1],
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
        rho_samples[ii,,tt+1] = metropolis_hastings_rho(alpha = alpha, n_items = n_items,
                                                                      rankings = aug_rankings[,,ii], metric = metric,
                                                                      rho = rho_samples[ii,,tt+1], leap_size = leap_size)

      }

      for (jj in 1:num_ranks){
        if(aug_method == "random"){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking(current_ranking = aug_rankings[jj,,ii],
                                                                               partial_ranking = R_obs[jj,,tt+1],
                                                                               alpha = alpha,
                                                                               rho = rho_samples[ii,,tt+1],
                                                                               n_items = n_items,
                                                                               metric = metric)
        }else if( (aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking_pseudo(current_ranking = aug_rankings[jj,,ii],
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
