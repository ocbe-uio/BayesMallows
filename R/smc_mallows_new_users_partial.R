#' @title SMC-Mallows new users partial
#'  @description Function to perform resample-move SMC algorithm where we receive new users with complete rankings
#' at each time step
#'   @param R_obs Matrix containing the full set of observed rankings of size n_assessors by n_items
#'   @param n_items Integer is the number of items in a ranking
#'   @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#'   @param leap_size leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution
#'   @param N Integer specifying the number of particles
#'   @param Time Integer specifying the number of time steps in the SMC algorithm
#'   @param logz_estimate Estimate of the partition function, computed with
#'   \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
#'   @param mcmc_kernel_app Interger value for the number of applications we apply the MCMC move kernel
#'   @param num_new_obs Integer value for the number of new observations (complete rankings) for each time step
#'   @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha
#'   @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the truncated exponential prior distribution of alpha.
#'   @param alpha_max  Maximum value of alpha in the truncated exponential
#'   prior distribution.
#'   @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
#' @return a set of particles each containing a value of rho and alpha
#' @export
smc_mallows_new_users_partial <- function(R_obs, n_items, metric, leap_size, N, Time, logz_estimate,
                                                 mcmc_kernel_app, num_new_obs, alpha_prop_sd,
                                                 lambda, alpha_max, aug_method){


  ######################
  ## Initialise Phase
  ######################
  n_users = dim(R_obs)[1] # this is total number of users

  # generate rho samples using uniform prior
  rho_samples = array(0,  c(N, n_items, (Time+1)))
  for (ii in 1:N){
    rho_samples[ii,,1] = sample(1:n_items, n_items, replace=FALSE)
  }

  # generate alpha samples using exponential prior
  alpha_samples = array(0, c(N, (Time+1)))
  alpha_samples[,1] = rexp(N, rate = 1)

  # this is to store the augmentations of the observed rankings for each particle
  aug_rankings = array(0, c(n_users, n_items, N))  # no. users by items by particles

  #########################
  ## New user situation
  #########################
  num_obs = 0

  for (tt in 1:Time){

    #print( paste("observe", tt, "out of" , Time))

    ###########################
    ## New Information
    ###########################
    # keep tally of how many ranking observations we have so far
    num_obs = num_obs + num_new_obs

    # create two ranking dataset to use for the reweight and move stages of the algorithm
    # Note:
    #new_observed_rankings = R_obs[((num_obs-num_new_obs+1):num_obs),]
    #all_observed_rankings = R_obs[(1:num_obs),]

    # propagate particles onto the next time step
    rho_samples[,,tt+1] = rho_samples[,,tt]
    alpha_samples[,tt+1] = alpha_samples[,tt]

    # calculate incremental weight and augmentation prob for each particle, based on new observed rankings
    log_inc_wgt = rep(0, N)

    ###########################
    # Augment partial rankings
    ###########################
    ranks = c(1:n_items)
    aug_prob = rep(1,N)

    for (ii in 1:N){
      for (jj in (num_obs-num_new_obs+1):num_obs){

        partial_ranking = R_obs[jj,]

        # find items missing from original observed ranking
        unranked_items = as.numeric(which(is.na(partial_ranking)))

        # find ranks missing from ranking
        missing_ranks = (ranks[!ranks %in% partial_ranking])

        # fill in missing ranks based on choice of augmentation method
        if(aug_method == "random"){

          # create new agumented ranking by sampling remaining ranks from set uniformly
          if(length(missing_ranks == 1)){partial_ranking[is.na(partial_ranking)] = missing_ranks}
          else {partial_ranking[is.na(partial_ranking)] <- sample(missing_ranks, size=length(missing_ranks), replace=F)}

          aug_rankings[jj,,ii] <- partial_ranking
          aug_prob[ii] = aug_prob[ii] * (1/factorial(length(missing_ranks)))

        }else if((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){

          # randomly permute the unranked items to give the order in which they will be allocated
          item_ordering = sample(unranked_items, size=length(unranked_items), replace=F)

          proposal = calculate_forward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                                   remaining_set = missing_ranks, rho = rho_samples[ii,,tt+1],
                                                   alpha = alpha_samples[ii,tt+1],
                                                   n_items = n_items, metric = metric)
          aug_rankings[jj,,ii] = proposal$aug_ranking
          aug_prob[ii] = aug_prob[ii] * proposal$forward_prob

        }else{

          print("Error: combined choice of metric and aug_method is incompatible")
          break()

        }
      }
    }

    ###########################
    # Re-weight
    ###########################
    for (ii in 1:N){
      # evaluate the log estimate of the partition function for a particular value of alpha
      log_z_alpha = BayesMallows:::get_partition_function(n_items = n_items, alpha =  alpha_samples[ii,tt+1],
                                                          logz_estimate = logz_estimate, metric = metric)
      log_likelihood = get_mallows_loglik(alpha = alpha_samples[ii,tt+1], rho = rho_samples[ii,,tt+1],
                                          n_items = n_items, rankings = aug_rankings[((num_obs-num_new_obs+1):num_obs),,ii],
                                          metric = metric)
      log_inc_wgt[ii] = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob[ii])
    }

    # normalise weights
    maxw = max(log_inc_wgt)
    w = exp(log_inc_wgt-maxw)
    norm_wgt = w/sum(w)


    ##############
    # Resample
    ##############
    # resample particles using multinomial resampling
    index = sample(1:N, prob=norm_wgt, size=N, rep=T)
    rho_samples[,,tt+1] = rho_samples[index,,tt+1]
    alpha_samples[,tt+1] = alpha_samples[index,tt+1]
    aug_rankings[(1:num_obs),,] = aug_rankings[(1:num_obs),,index]

    ##################
    # Move step
    ##################
    for (ii in 1:N){
      for (kk in 1:mcmc_kernel_app){
        # move each particle containing sample of rho and alpha by using the MCMC kernels
        rho_samples[ii,,tt+1] = metropolis_hastings_rho(alpha = alpha_samples[ii,tt+1], n_items = n_items,
                                                        rankings = aug_rankings[(1:num_obs),,ii],
                                                        metric = metric, rho = rho_samples[ii,,tt+1],
                                                        leap_size = leap_size)

      }
      alpha_samples[ii,tt+1] = metropolis_hastings_alpha_update(alpha = alpha_samples[ii,tt+1], n_items = n_items,
                                                         rankings = aug_rankings[(1:num_obs),,ii],
                                                         metric = metric, rho = rho_samples[ii,,tt+1],
                                                         logz_estimate = logz_estimate, alpha_prop_sd = alpha_prop_sd,
                                                         lambda = lambda, alpha_max = alpha_max)
      for (jj in 1:num_obs){
        if (aug_method == "random"){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking(R_curr = aug_rankings[jj,,ii],
                                                                 R_obs = R_obs[jj,], alpha = alpha_samples[ii,tt+1],
                                                                 rho = rho_samples[ii,,tt+1], n_items = n_items,
                                                                 metric = metric)

        }else if((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking_pseudo(alpha = alpha_samples[ii,tt+1], rho = rho_samples[ii,,tt+1],
                                                                      n_items = n_items, partial_ranking = R_obs[jj,],
                                                                      current_ranking = aug_rankings[jj,,ii], metric = metric)
        }
      }
    }

  }
  # return the history of the particles and their values
  smc_list <- list("rho_samples" = rho_samples, "alpha_samples" = alpha_samples)
  return(smc_list)
}






#' @title SMC-mallows new users partial (alpha fixed)
#' @description Function to perform resample-move SMC algorithm where we receive new users with complete rankings
#' at each time step
#'   @param R_obs Matrix containing the full set of observed rankings of size n_assessors by n_items
#'   @param n_items Integer is the number of items in a ranking
#'   @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#'   @param leap_size leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution
#'   @param N Integer specifying the number of particles
#'   @param Time Integer specifying the number of time steps in the SMC algorithm
#'   @param logz_estimate Estimate of the partition function, computed with
#'   \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
#'   @param mcmc_kernel_app Interger value for the number of applications we apply the MCMC move kernel
#'   @param num_new_obs Integer value for the number of new observations (complete rankings) for each time step
#'   @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"
#'   @param alpha A numeric value of the scale parameter which is known and fixed
#' @return a set of particles each containing a value of rho and alpha
#' @export
smc_mallows_new_users_partial_alpha_fixed <- function(R_obs, n_items, metric, leap_size, N, Time, logz_estimate,
                                          mcmc_kernel_app, num_new_obs, aug_method, alpha){


  ######################
  ## Initialise Phase
  ######################
  n_users = dim(R_obs)[1] # this is total number of users

  # generate rho samples using uniform prior
  rho_samples = array(0,  c(N, n_items, (Time+1)))
  for (ii in 1:N){
    rho_samples[ii,,1] = sample(1:n_items, n_items, replace=FALSE)
  }

  # this is to store the augmentations of the observed rankings for each particle
  aug_rankings = array(0, c(n_users, n_items, N))  # no. users by items by particles

  #########################
  ## New user situation
  #########################
  num_obs = 0

  for (tt in 1:Time){

    #print( paste("observe", tt, "out of" , Time))

    ###########################
    ## New Information
    ###########################
    # keep tally of how many ranking observations we have so far
    num_obs = num_obs + num_new_obs

    # create two ranking dataset to use for the reweight and move stages of the algorithm
    # Note:
    #new_observed_rankings = R_obs[((num_obs-num_new_obs+1):num_obs),]
    #all_observed_rankings = R_obs[(1:num_obs),]

    # propagate particles onto the next time step
    rho_samples[,,tt+1] = rho_samples[,,tt]

    # calculate incremental weight and augmentation prob for each particle, based on new observed rankings
    log_inc_wgt = rep(0, N)

    ###########################
    # Augment partial rankings
    ###########################
    ranks = c(1:n_items)
    aug_prob = rep(1,N)

    for (ii in 1:N){
      for (jj in (num_obs-num_new_obs+1):num_obs){

        partial_ranking = R_obs[jj,]

        # find items missing from original observed ranking
        unranked_items = as.numeric(which(is.na(partial_ranking)))

        # find ranks missing from ranking
        missing_ranks = (ranks[!ranks %in% partial_ranking])

        # fill in missing ranks based on choice of augmentation method
        if(aug_method == "random"){

          # create new agumented ranking by sampling remaining ranks from set uniformly
          if(length(missing_ranks == 1)){partial_ranking[is.na(partial_ranking)] = missing_ranks}
          else {partial_ranking[is.na(partial_ranking)] <- sample(missing_ranks, size=length(missing_ranks), replace=F)}

          aug_rankings[jj,,ii] <- partial_ranking
          aug_prob[ii] = aug_prob[ii] * (1/factorial(length(missing_ranks)))

        }else if((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){

          # randomly permute the unranked items to give the order in which they will be allocated
          item_ordering = sample(unranked_items, size=length(unranked_items), replace=F)

          proposal = calculate_forward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                                   remaining_set = missing_ranks, rho = rho_samples[ii,,tt+1],
                                                   alpha = alpha, n_items = n_items, metric = metric)
          aug_rankings[jj,,ii] = proposal$aug_ranking
          aug_prob[ii] = aug_prob[ii] * proposal$forward_prob

        }else{

          print("Error: combined choice of metric and aug_method is incompatible")
          break()

        }
      }
    }

    ###########################
    # Re-weight
    ###########################
    for (ii in 1:N){
      # evaluate the log estimate of the partition function for a particular value of alpha
      log_z_alpha = BayesMallows:::get_partition_function(n_items = n_items, alpha =  alpha,
                                                          logz_estimate = logz_estimate, metric = metric)
      log_likelihood = get_mallows_loglik(alpha = alpha, rho = rho_samples[ii,,tt+1],
                                          n_items = n_items, rankings = aug_rankings[((num_obs-num_new_obs+1):num_obs),,ii],
                                          metric = metric)
      log_inc_wgt[ii] = log_likelihood - num_new_obs*log_z_alpha - log(aug_prob[ii])
    }

    # normalise weights
    maxw = max(log_inc_wgt)
    w = exp(log_inc_wgt-maxw)
    norm_wgt = w/sum(w)


    ##############
    # Resample
    ##############
    # resample particles using multinomial resampling
    index = sample(1:N, prob=norm_wgt, size=N, rep=T)
    rho_samples[,,tt+1] = rho_samples[index,,tt+1]
    aug_rankings[(1:num_obs),,] = aug_rankings[(1:num_obs),,index]

    ##################
    # Move step
    ##################
    for (ii in 1:N){
      for (kk in 1:mcmc_kernel_app){
        # move each particle containing sample of rho and alpha by using the MCMC kernels
        rho_samples[ii,,tt+1] = metropolis_hastings_rho(alpha = alpha, n_items = n_items,
                                                        rankings = aug_rankings[(1:num_obs),,ii],
                                                        metric = metric, rho = rho_samples[ii,,tt+1],
                                                        leap_size = leap_size)

      }

      for (jj in 1:num_obs){
        if (aug_method == "random"){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking(R_curr = aug_rankings[jj,,ii],
                                                                 R_obs = R_obs[jj,], alpha = alpha,
                                                                 rho = rho_samples[ii,,tt+1], n_items = n_items,
                                                                 metric = metric)

        }else if((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman")) ){
          aug_rankings[jj,,ii] = metropolis_hastings_aug_ranking_pseudo(alpha = alpha, rho = rho_samples[ii,,tt+1],
                                                                        n_items = n_items, partial_ranking = R_obs[jj,],
                                                                        current_ranking = aug_rankings[jj,,ii], metric = metric)
        }
      }
    }

  }
  # return the history of the particles and their values
  smc_list <- list("rho_samples" = rho_samples)
  return(smc_list)
}
