##############################################################################
#SMC-Mallows using New Users with Complete Rankings as Observed Data
##############################################################################

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
source("post_processing_functions.R")


#################################################################################################
############################### SCRIPT ##########################################################
#################################################################################################

smc_mallows_new_users_complete<- function(R_obs, n_items, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app, num_new_obs){
  
  # @description Function to perform resample-move SMC algorithm where we receive new users with complete rankings at each time step
  
  # INPUT:
  #   @param R_obs Matrix containing the full set of observed rankings of size n_assessors by n_items
  #   @param n_items Integer is the number of items in a ranking
  #   @param metric A character string specifying the distance metric to use in the
  #   Bayesian Mallows Model. Available options are \code{"footrule"},
  #   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
  #   \code{"ulam"}.
  #   @param leap_size leap_size Integer specifying the step size of the leap-and-shift
  #   proposal distribution
  #   @param N Integer specifying the number of particles
  #   @param Time Integer specifying the number of time steps in the SMC algorithm
  #   @param logz_estimate Estimate of the partition function, computed with
  #   \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
  #   @param mcmc_kernel_app Interger value for the number of applications we apply the MCMC move kernel
  #   @param num_new_obs Integer value for the number of new observations (complete rankings) for each time step
  
  # OUTPUT: a set of particles each containing a value of rho and alpha
  
  ######################
  ## Initialise Phase
  ######################
  n_users = dim(R_obs)[1] # this is total number of users
  Time = dim(R_obs)[1]/num_new_obs # determine number of time steps for for loop
  
  # generate rho samples using uniform prior
  rho_samples = array(0,  c(N, n_items, (n_users+Time+1)))
  for (ii in 1:N){
    rho_samples[ii,,1] = sample(1:n_items, n_items, replace=FALSE)
  }
  
  # generate alpha samples using exponential prior
  alpha_samples = array(0, c(N, (n_users+Time+1)))
  alpha_samples[,1] = rexp(N, rate = 1)
  
  #########################
  ## New user situation
  #########################
  num_obs = 0
  
  for (tt in 1:Time){
    print( paste("observe", tt, "out of" , Time)) 
    
    # keep tally of how many ranking observations we have so far
    num_obs = num_obs + num_new_obs
    
    ###########################
    ## New Information
    ###########################
    # create two rannking dataset to use for the reweight and move stages of the algorithm
    new_observed_rankings = R_obs[((num_obs-num_new_obs+1):num_obs),]
    all_observed_rankings = R_obs[(1:num_obs),]
    
    # propagate particles onto the next time step
    rho_samples[,,tt+1] = rho_samples[,,tt]
    alpha_samples[,tt+1] = alpha_samples[,tt]
    
    ###########################
    # Re-weight
    ###########################
    # calculate incremental weight for each particle, based on new observed rankings
    log_inc_wgt = rep(0, N) 
    
    for (ii in 1:N){ 
        # evaluate the log estimate of the partition function for a particular value of alpha
        log_z_alpha = BayesMallows:::get_partition_function(n_items = n_items, alpha =  alpha_samples[ii,tt+1], 
                                                            logz_estimate = logz_estimate, metric = metric)
        log_likelihood = get_mallows_loglik(alpha = alpha_samples[ii,tt+1], rho = rho_samples[ii,,tt+1], 
                                            n_items = n_items, rankings = new_observed_rankings, metric = metric) 
        log_inc_wgt[ii] = log_likelihood - num_new_obs*log_z_alpha 
    }
    
    # normalise weights
    maxw = max(log_inc_wgt)
    w = exp(log_inc_wgt-maxw)
    norm_wgt = w/sum(w) 
    
    ##############
    # Resample
    ##############
    # resample particles using multinomial resampling 
    index = sample(1:N,prob=norm_wgt,size=N,rep=T)   
    rho_samples[,,tt+1] = rho_samples[index,,tt+1]
    alpha_samples[,tt+1] = alpha_samples[index,tt+1]

    ################## 
    # Move step
    ##################
    for (ii in 1:N){
      for (kk in 1:mcmc_kernel_app){
        # move each particle containing sample of rho and alpha by using the MCMC kernels
        rho_samples[ii,,tt+1] = metropolis_hastings_rho(alpha = alpha_samples[ii,tt+1], n_items = n_items, 
                                                        rankings = all_observed_rankings, 
                                                        metric = metric, rho = rho_samples[ii,,tt+1], 
                                                        leap_size = leap_size)
        alpha_samples[ii,tt+1] = metropolis_hastings_alpha(alpha = alpha_samples[ii,tt+1], n_items = n_items, 
                                                           rankings = all_observed_rankings, 
                                                           metric = metric, rho = rho_samples[ii,,tt+1], 
                                                           logz_estimate = logz_estimate)

      }
      

    }
    
  }
  # return the history of the particles and their values 
  smc_list <- list("rho_samples" = rho_samples, "alpha_samples" = alpha_samples)  
  return(smc_list)
}


