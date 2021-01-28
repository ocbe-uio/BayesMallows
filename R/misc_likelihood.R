#/' Log-likelihood evaluation for a Mallows model
#/'
#/' @description Compute the log-likelihood value of the Mallows model parameters for a ranking dataset.
#/' @param rho Numeric vector of length n_items. \code{rho} must be a permutation of the first n_items integers and corresponds to the location of the Mallows distribution. The location coincides with the mode (most probable permutation).
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Kendall distance.
#/' @param metric Character string specifying the distance measure to use. Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"}, \code{"ulam"} for \code{n_items<=95}, \code{"footrule"} for \code{n_items<=50} and \code{"spearman"} for \code{n_items<=14}.
#/' @param rankings A matrix with observed rankings in each row.
#/' @param obs_freq A vector of observation frequencies (weights) to apply to each row in \code{rankings}.
#/'   This can speed up computation if a large number of assessors share the same
#/'   rank pattern.
#/' @return The log-likelihood value corresponding to one or more observed rankings under the Mallows rank model with distance specified by the \code{metric} argument.
#/'

log_lik_db <- function(rho,alpha,metric,rankings,obs_freq){

  N <- sum(obs_freq)
  n_items <- ncol(rankings)

  if(metric%in%c("kendall","cayley","hamming")){
    log_lik <- -(alpha*rank_dist_sum(rankings=t(rankings),rho=rho,metric=metric,obs_freq=obs_freq)+N*get_partition_function(n_items=n_items,alpha=alpha*n_items,metric=metric))
  }

  if(metric%in%c("ulam","footrule","spearman")){
    pfd <- dplyr::filter(partition_function_data,
                         .data$metric == !!metric, .data$n_items == !!n_items,
                         .data$type == "cardinalities")
    if(nrow(pfd) == 0){
      stop("Given number of items currently not available for the specified metric")
    } else{
      card <- pfd$values[[1]]
    }
    log_lik <- -(alpha*rank_dist_sum(rankings=t(rankings),rho=rho,metric=metric,obs_freq=obs_freq)+N*get_partition_function(alpha=alpha*n_items,n_items=n_items,metric=metric,cardinalities=card))
  }

  return(log_lik)
}



#/' Log-likelihood evaluation for a Mallows mixture model
#/'
#/' @description Compute the log-likelihood value of the Mallows mixture model parameters for a ranking dataset.
#/' @param rho A matrix of size \code{n_clusters x n_items} whose rows are permutations of the first n_items integers corresponding to the modal rankings of the Mallows mixture components.
#/' @param alpha A vector of \code{n_clusters} non-negative scalars specifying the scale (precision) parameters of the Mallows mixture components.
#/' @param weights A vector of \code{n_clusters} non-negative scalars specifying the mixture weights.
#/' @param metric Character string specifying the distance measure to use.
#/' @param rankings A matrix with observed rankings in each row.
#/' @param obs_freq A vector of observation frequencies (weights) to apply to each row in \code{rankings}.
#/'   This can speed up computation if a large number of assessors share the same
#/'   rank pattern. If \code{NULL}, it means that each row of
#/'   \code{rankings} is multiplied by 1. If provided, \code{obs_freq} must have
#/'   the same number of elements as there are rows in \code{rankings}, and
#/'   \code{rankings} cannot be \code{NULL}.
#/'
#/' @return The log-likelihood value corresponding to one or more observed rankings under the Mallows mixture model with distance specified by the \code{metric} argument.
#/'

log_lik_db_mix <- function(rho,alpha,weights,metric,rankings,obs_freq){

  L <- length(obs_freq)
  n_clusters <- length(weights)
  temp <- matrix(NA,nrow=n_clusters,ncol=L)
  for(l in 1:L){
    for(g in 1:n_clusters){
      temp[g,l] <- exp(log_lik_db(rho=rho[g,],alpha=alpha[g],metric=metric,rankings=rankings[l,,drop=FALSE],obs_freq=obs_freq[l]))
    }
  }
  log_lik <- sum(log(weights%*%temp))
  return(log_lik)
}

