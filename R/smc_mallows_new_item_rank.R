smc_mallows_new_item_rank <- function(n_items, R_obs, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app,
                                      alpha_prop_sd, lambda, alpha_max, aug_method) {

  # @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user
  # at each time step. Each correction and augmentation is done by filling in the missing item ranks using pseudlikelihood augmentation.

  # INPUT:
  #   @param n_items Integer is the number of items in a ranking
  #   @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time time steps
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
  #   @param alpha_prop_sd Numeric value of the standard deviation of the prior distribution for alpha
  #   @param lambda Strictly positive numeric value specifying the rate parameter
  #   of the truncated exponential prior distribution of alpha.
  #   @param alpha_max  Maximum value of alpha in the truncated exponential
  #   prior distribution.
  #   @param aug_method A character string specifying the approach for filling in the missing data, options are "pseudolikelihood" or "random"


  # OUTPUT: a 3d matrix containing the samples of rho and alpha from the SMC algorithm

  ######################
  ## Initialise Phase
  ######################
  # Generate N initial samples of rho using the uniform prior
  rho_samples <- array(0, c(N, n_items, Time + 1))
  for (ii in 1:N) {
    rho_samples[ii, , 1] <- sample(1:n_items, n_items, replace = FALSE)
  }

  alpha_samples <- array(0, c(N, Time + 1))
  alpha_samples[, 1] <- rexp(N, rate = 1)

  ######################
  ## Augment Rankings
  ######################
  num_ranks <- dim(R_obs[, , 1])[1]

  # each particle has its own set of augmented rankings
  aug_rankings <- array(0, c(num_ranks, n_items, N))
  prev_aug_rankings <- array(0, c(num_ranks, n_items, N))


  # augment incomplete ranks to initialise
  ranks <- c(1:n_items)

  # total correction prob
  total_correction_prob <- rep(1, N)


  # iterate through each observed ranking and create new "corrected" augmented rankings
  for (ii in 1:N) {
    # set t-1 generation to old as we sample for t new
    prev_aug_rankings[, , ii] <- aug_rankings[, , ii]

    # make the correction
    for (jj in 1:num_ranks) {


      # fill in missing ranks based on choice of augmentation method
      if (aug_method == "random") {

        # find elements missing from original observed ranking
        partial_ranking <- R_obs[jj, , 1]

        remaining_set <- (ranks[!ranks %in% partial_ranking])

        # create new agumented ranking by sampling remaining ranks from set uniformly
        if (length(remaining_set == 1)) {
          partial_ranking[is.na(partial_ranking)] <- remaining_set
        }
        else {
          partial_ranking[is.na(partial_ranking)] <- sample(remaining_set, size = length(remaining_set), replace = F)
        }

        aug_rankings[jj, , ii] <- partial_ranking
        total_correction_prob[ii] <- total_correction_prob[ii] * (1 / factorial(length(remaining_set)))
      } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {


        # find items missing from original observed ranking
        unranked_items <- as.numeric(which(is.na(R_obs[jj, , 1])))

        # find unallocated ranks from original observed ranking
        remaining_set <- ranks[!ranks %in% R_obs[jj, , 1]]

        # randomly permute the unranked items to give the order in which they will be allocated
        item_ordering <- sample(unranked_items, size = length(unranked_items), replace = F)

        proposal <- calculate_forward_probability(
          item_ordering = item_ordering, partial_ranking = R_obs[jj, , 1],
          remaining_set = remaining_set, rho = rho_samples[ii, , 1], alpha = alpha_samples[ii, 1],
          n_items = n_items, metric = metric
        )
        aug_rankings[jj, , ii] <- proposal$aug_ranking
        total_correction_prob[ii] <- total_correction_prob[ii] * proposal$forward_prob
      } else {
        stop(
          "Combined choice of metric and aug_method is incompatible. ",
          "The value is TRUE, so the script must end here"
        )
      }
    }
  }

  ###########################
  # Re-weight
  ###########################
  # incremental weight for each particle, based on new observed rankings
  log_inc_wgt <- rep(0, N)

  for (ii in 1:N) {
    # evaluate the log estimate of the partition function for a particular value of alpha
    log_z_alpha <- get_partition_function(
      n_items = n_items, alpha = alpha_samples[ii, 1],
      logz_estimate = logz_estimate, metric = metric
    )
    log_likelihood <- get_mallows_loglik(
      alpha = alpha_samples[ii, 1], rho = rho_samples[ii, , 1],
      n_items = n_items, rankings = aug_rankings[, , ii], metric = metric
    )
    log_inc_wgt[ii] <- log_likelihood - num_ranks * log_z_alpha - log(total_correction_prob[ii])
  }

  # update weights
  maxw <- max(log_inc_wgt)
  w <- exp(log_inc_wgt - maxw)
  norm_wgt <- w / sum(w)


  ##############
  # Resample
  ##############
  index <- sample(1:N, prob = norm_wgt, size = N, replace = T)
  rho_samples[, , 1] <- rho_samples[index, , 1]
  alpha_samples[, 1] <- alpha_samples[index, 1]
  aug_rankings[, , ] <- aug_rankings[, , index]


  ##################
  # Move step
  ##################
  for (ii in 1:N) {
    for (kk in 1:mcmc_kernel_app) {
      rho_samples[ii, , 1] <- metropolis_hastings_rho(
        alpha = alpha_samples[ii, 1], n_items = n_items,
        rankings = aug_rankings[, , ii], metric = metric,
        rho = rho_samples[ii, , 1], leap_size = leap_size
      )

      alpha_samples[ii, 1] <- metropolis_hastings_alpha(
        alpha = alpha_samples[ii, 1], n_items = n_items,
        rankings = aug_rankings[, , ii], metric = metric,
        rho = rho_samples[ii, , 1],
        logz_estimate = logz_estimate, alpha_prop_sd = alpha_prop_sd,
        lambda = lambda, alpha_max = alpha_max
      )
    }
    for (jj in 1:num_ranks) {
      if (aug_method == "random") {
        aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking(
          alpha = alpha_samples[ii, 1],
          rho = rho_samples[ii, , 1], n_items = n_items,
          R_obs = R_obs[jj, , 1],
          R_curr = aug_rankings[jj, , ii], metric = metric
        )
      } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
        aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking_pseudo(
          alpha = alpha_samples[ii, 1],
          rho = rho_samples[ii, , 1], n_items = n_items,
          partial_ranking = R_obs[jj, , 1],
          current_ranking = aug_rankings[jj, , ii], metric = metric
        )
      }
    }
  }

  #########################
  ## Loop for t=1,...,Time
  #########################
  for (tt in 1:(Time - 1)) {
    message("iteration ", tt, " out of ", Time)

    ###########################
    ## New Information
    ###########################
    # new observed item ranks from each user, need to update augmented rankings
    rho_samples[, , tt + 1] <- rho_samples[, , tt]
    alpha_samples[, tt + 1] <- alpha_samples[, tt]

    # total correction prob
    particle_correction_prob <- rep(1, N)

    # iterate through each observed ranking and create new "corrected" augmented rankings
    for (ii in 1:N) {
      # set t-1 generation to old as we sample for t new
      prev_aug_rankings[, , ii] <- aug_rankings[, , ii]

      # make the correction
      for (jj in 1:num_ranks) {
        if (aug_method == "random") {
          check_correction <- correction_kernel(
            R_curr = aug_rankings[jj, , ii],
            R_obs = R_obs[jj, , tt + 1], n_items = n_items
          )
          aug_rankings[jj, , ii] <- check_correction$ranking

          particle_correction_prob[ii] <- particle_correction_prob[ii] * check_correction$correction_prob
        } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
          check_correction <- correction_kernel_pseudo(
            R_curr = aug_rankings[jj, , ii],
            R_obs = R_obs[jj, , tt + 1], rho = rho_samples[ii, , tt + 1],
            alpha = alpha_samples[ii, tt + 1], n_items = n_items, metric = metric
          )

          aug_rankings[jj, , ii] <- check_correction$ranking

          # these probs are in real scale
          particle_correction_prob[ii] <- particle_correction_prob[ii] * check_correction$correction_prob
        } else {
          stop("Combined choice of metric and aug_method is incompatible")
        }
      }
    }



    ###########################
    # Re-weight
    ###########################
    # incremental weight for each particle, based on new observed rankings
    log_inc_wgt <- rep(0, N)

    for (ii in 1:N) {
      log_inc_wgt[ii] <- get_mallows_loglik(
        alpha = alpha_samples[ii, tt + 1], rho = rho_samples[ii, , tt + 1],
        n_items = n_items, rankings = aug_rankings[, , ii], metric = metric
      ) -
        get_mallows_loglik(
          alpha = alpha_samples[ii, tt + 1], rho = rho_samples[ii, , tt + 1],
          n_items = n_items, rankings = prev_aug_rankings[, , ii], metric = metric
        ) -
        log(particle_correction_prob[ii])
    }

    # update unnormalised weights
    maxw <- max(log_inc_wgt)
    w <- exp(log_inc_wgt - maxw)
    norm_wgt <- w / sum(w)

    ##############
    # Resample
    ##############
    index <- sample(1:N, prob = norm_wgt, size = N, replace = T)
    rho_samples[, , tt + 1] <- rho_samples[index, , tt + 1]
    alpha_samples[, tt + 1] <- alpha_samples[index, tt + 1]
    aug_rankings[, , ] <- aug_rankings[, , index]

    ##################
    # Move step
    ##################
    for (ii in 1:N) {
      for (kk in 1:mcmc_kernel_app) {
        rho_samples[ii, , tt + 1] <- metropolis_hastings_rho(
          alpha = alpha_samples[ii, tt + 1], n_items = n_items,
          rankings = aug_rankings[, , ii], metric = metric,
          rho = rho_samples[ii, , tt + 1], leap_size = leap_size
        )

        alpha_samples[ii, tt + 1] <- metropolis_hastings_alpha(
          alpha = alpha_samples[ii, tt + 1], n_items = n_items,
          rankings = aug_rankings[, , ii], metric = metric,
          rho = rho_samples[ii, , tt + 1],
          logz_estimate = logz_estimate, alpha_prop_sd = alpha_prop_sd,
          lambda = lambda, alpha_max = alpha_max
        )
      }
      for (jj in 1:num_ranks) {
        if (aug_method == "random") {
          aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking(
            R_curr = aug_rankings[jj, , ii],
            R_obs = R_obs[jj, , tt + 1],
            alpha = alpha_samples[ii, tt + 1],
            rho = rho_samples[ii, , tt + 1], n_items = n_items,
            metric = metric
          )
        } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
          aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking_pseudo(
            current_ranking = aug_rankings[jj, , ii],
            partial_ranking = R_obs[jj, , tt + 1],
            alpha = alpha_samples[ii, tt + 1],
            rho = rho_samples[ii, , tt + 1], n_items = n_items,
            metric = metric
          )
        }
      }
    }
  }

  ############################
  # Post Processing
  ############################
  smc_list <- list("rho_samples" = rho_samples, "alpha_samples" = alpha_samples)
  return(smc_list)
}

smc_mallows_new_item_rank_alpha_fixed <- function(alpha, n_items, R_obs, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app,
                                                  alpha_prop_sd, lambda, alpha_max, aug_method) {

  # @description Function to perform resample-move SMC algorithm where we receive a new item ranks from an existing user
  # at each time step. Each correction and augmentation is done by filling in the missing item ranks randomly.

  # INPUT:
  #   @param alpha A numeric value of the true scale parameter
  #   @param n_items Integer is the number of items in a ranking
  #   @param R_obs 3D matrix of size n_assessors by n_items by Time containing a set of observed rankings of Time time steps
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
  #   @param mcmc_kernel_app Integer value for the number of applications we apply the MCMC move kernel

  # OUTPUT: a 3d matrix containing the samples of rho and alpha from the SMC algorithm

  ######################
  ## Initialise Phase
  ######################
  # Generate N initial samples of rho using the uniform prior
  rho_samples <- array(0, c(N, n_items, Time + 1))
  for (ii in 1:N) {
    rho_samples[ii, , 1] <- sample(1:n_items, n_items, replace = FALSE)
  }

  ######################
  ## Augment Rankings
  ######################
  num_ranks <- dim(R_obs[, , 1])[1]

  # each particle has its own set of augmented rankings
  aug_rankings <- array(0, c(num_ranks, n_items, N))
  prev_aug_rankings <- array(0, c(num_ranks, n_items, N))

  # augment incomplete ranks to initialise
  ranks <- c(1:n_items)

  # total correction prob
  total_correction_prob <- rep(1, N)


  # iterate through each observed ranking and create new "corrected" augmented rankings
  for (ii in 1:N) {
    # set t-1 generation to old as we sample for t new
    prev_aug_rankings[, , ii] <- aug_rankings[, , ii]

    # make the correction
    for (jj in 1:num_ranks) {


      # fill in missing ranks based on choice of augmentation method
      if (aug_method == "random") {

        # find elements missing from original observed ranking
        partial_ranking <- R_obs[jj, , 1]

        remaining_set <- (ranks[!ranks %in% partial_ranking])

        # create new agumented ranking by sampling remaining ranks from set uniformly
        if (length(remaining_set == 1)) {
          partial_ranking[is.na(partial_ranking)] <- remaining_set
        }
        else {
          partial_ranking[is.na(partial_ranking)] <- sample(remaining_set, size = length(remaining_set), replace = F)
        }

        aug_rankings[jj, , ii] <- partial_ranking
        total_correction_prob[ii] <- total_correction_prob[ii] * (1 / factorial(length(remaining_set)))
      } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {


        # find items missing from original observed ranking
        unranked_items <- as.numeric(which(is.na(R_obs[jj, , 1])))

        # find unallocated ranks from original observed ranking
        remaining_set <- ranks[!ranks %in% R_obs[jj, , 1]]

        # randomly permute the unranked items to give the order in which they will be allocated
        item_ordering <- sample(unranked_items, size = length(unranked_items), replace = F)

        proposal <- calculate_forward_probability(
          item_ordering = item_ordering, partial_ranking = R_obs[jj, , 1],
          remaining_set = remaining_set, rho = rho_samples[ii, , 1], alpha = alpha,
          n_items = n_items, metric = metric
        )
        aug_rankings[jj, , ii] <- proposal$aug_ranking
        total_correction_prob[ii] <- total_correction_prob[ii] * proposal$forward_prob
      } else {
        stop("Combined choice of metric and aug_method is incompatible")
      }
    }
  }

  ###########################
  # Re-weight
  ###########################
  log_inc_wgt <- rep(0, N)

  for (ii in 1:N) {
    # evaluate the log estimate of the partition function for a particular value of alpha
    log_z_alpha <- get_partition_function(
      n_items = n_items, alpha = alpha,
      logz_estimate = logz_estimate, metric = metric
    )
    log_likelihood <- get_mallows_loglik(
      alpha = alpha, rho = rho_samples[ii, , 1],
      n_items = n_items, rankings = aug_rankings[, , ii], metric = metric
    )
    log_inc_wgt[ii] <- log_likelihood - num_ranks * log_z_alpha - log(total_correction_prob[ii])
  }

  # update weights
  maxw <- max(log_inc_wgt)
  w <- exp(log_inc_wgt - maxw)
  norm_wgt <- w / sum(w)

  ##############
  # Resample
  ##############
  index <- sample(1:N, prob = norm_wgt, size = N, replace = T)
  rho_samples[, , 1] <- rho_samples[index, , 1]
  aug_rankings[, , ] <- aug_rankings[, , index]

  ##################
  # Move step
  ##################
  for (ii in 1:N) {
    for (kk in 1:mcmc_kernel_app) {
      rho_samples[ii, , 1] <- metropolis_hastings_rho(
        alpha = alpha, n_items = n_items,
        rankings = aug_rankings[, , ii], metric = metric,
        rho = rho_samples[ii, , 1], leap_size = leap_size
      )
    }
    for (jj in 1:num_ranks) {
      if (aug_method == "random") {
        aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking(
          alpha = alpha,
          rho = rho_samples[ii, , 1], n_items = n_items,
          R_obs = R_obs[jj, , 1],
          R_curr = aug_rankings[jj, , ii], metric = metric
        )
      } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
        aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking_pseudo(
          alpha = alpha,
          rho = rho_samples[ii, , 1], n_items = n_items,
          partial_ranking = R_obs[jj, , 1],
          current_ranking = aug_rankings[jj, , ii], metric = metric
        )
      }
    }
  }

  #########################
  ## Loop for t=1,...,Time
  #########################

  for (tt in 1:(Time - 1)) {
    message("We are now on iteration ", tt, " out of ", Time)

    ###########################
    ## New Information
    ###########################
    # new observed item ranks from each user, need to update augmented rankings
    rho_samples[, , tt + 1] <- rho_samples[, , tt]

    # total correction prob
    total_correction_prob <- rep(1, N)

    # iterate through each observed ranking and create new "corrected" augmented rankings
    for (ii in 1:N) {
      # set t-1 generation to old as we sample for t new
      prev_aug_rankings[, , ii] <- aug_rankings[, , ii]

      # make the correction
      for (jj in 1:num_ranks) {
        if (aug_method == "random") {
          check_correction <- correction_kernel(
            R_curr = aug_rankings[jj, , ii],
            R_obs = R_obs[jj, , tt + 1], n_items = n_items
          )

          aug_rankings[jj, , ii] <- check_correction$ranking

          total_correction_prob[ii] <- total_correction_prob[ii] * check_correction$correction_prob
        } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
          check_correction <- correction_kernel_pseudo(
            R_curr = aug_rankings[jj, , ii],
            R_obs = R_obs[jj, , tt + 1], rho = rho_samples[ii, , tt + 1],
            alpha = alpha, n_items = n_items, metric = metric
          )

          aug_rankings[jj, , ii] <- check_correction$ranking

          # these probs are in real scale
          total_correction_prob[ii] <- total_correction_prob[ii] * check_correction$correction_prob
        } else {
          stop("Combined choice of metric and aug_method is incompatible")
        }
      }
    }

    ###########################
    # Re-weight
    ###########################
    # incremental weight for each particle, based on new observed rankings
    log_inc_wgt <- rep(0, N)

    for (ii in 1:N) {
      log_inc_wgt[ii] <- get_mallows_loglik(
        alpha = alpha, rho = rho_samples[ii, , tt + 1],
        n_items = n_items, rankings = aug_rankings[, , ii], metric = metric
      ) -
        get_mallows_loglik(
          alpha = alpha, rho = rho_samples[ii, , tt + 1],
          n_items = n_items, rankings = prev_aug_rankings[, , ii], metric = metric
        ) -
        log(total_correction_prob[ii])
    }

    # update unnormalised weights
    maxw <- max(log_inc_wgt)
    w <- exp(log_inc_wgt - maxw)
    norm_wgt <- w / sum(w)

    ##############
    # Resample
    ##############
    index <- sample(1:N, prob = norm_wgt, size = N, replace = T)
    rho_samples[, , tt + 1] <- rho_samples[index, , tt + 1]

    ##################
    # Move step
    ##################
    for (ii in 1:N) {
      for (kk in 1:mcmc_kernel_app) {
        rho_samples[ii, , tt + 1] <- metropolis_hastings_rho(
          alpha = alpha, n_items = n_items,
          rankings = aug_rankings[, , ii], metric = metric,
          rho = rho_samples[ii, , tt + 1], leap_size = leap_size
        )
      }
      for (jj in 1:num_ranks) {
        if (aug_method == "random") {
          aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking(
            alpha = alpha,
            rho = rho_samples[ii, , tt + 1], n_items = n_items,
            R_obs = R_obs[jj, , tt + 1],
            R_curr = aug_rankings[jj, , ii], metric = metric
          )
        } else if ((aug_method == "pseudolikelihood") & (metric %in% c("footrule", "spearman"))) {
          aug_rankings[jj, , ii] <- metropolis_hastings_aug_ranking_pseudo(
            alpha = alpha,
            rho = rho_samples[ii, , tt + 1], n_items = n_items,
            partial_ranking = R_obs[jj, , tt + 1],
            current_ranking = aug_rankings[jj, , ii], metric = metric
          )
        }
      }
    }
  }
  ############################
  # Post Processing
  ############################
  smc_list <- list("rho_samples" = rho_samples)
  return(smc_list)
}
