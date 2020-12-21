#/' exp_d_tau
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Kendall distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Kendall metric under the Mallows rank model with the Kendall distance.

exp_d_tau <- function(alpha,n_items){
  if(alpha>0){
    idx <- 1:n_items
    out <- n_items*exp(-alpha)/(1-exp(-alpha))-sum((idx*exp(-idx*alpha))/(1-exp(-idx*alpha)))
  }else{
    if(alpha==0){
      out <- n_items*(n_items-1)/4
    }else{
      stop("alpha must be a non-negative value")
    }
  }
  return(out)
}

#/' exp_d_cay
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Cayley distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Cayley metric under the Mallows rank model with the Cayley distance.

exp_d_cay <- function(alpha,n_items){
  idx <- 1:(n_items-1)
  out <- sum(idx/(idx+exp(alpha)))
  return(out)
}

#/' exp_d_ham
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Hamming distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Hamming metric under the Mallows rank model with the Hamming distance.

exp_d_ham <- function(alpha,n_items){
  idx <- 0:n_items
  out <- n_items-exp(alpha)*sum(((exp(alpha)-1)^idx[-(n_items+1)])/base::factorial(idx[-(n_items+1)]))/sum(((exp(alpha)-1)^idx)/base::factorial(idx))
  return(out)
}

#/' exp_d_ulam
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Ulam distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Ulam metric under the Mallows rank model with the Ulam distance.

# The function based on the command from the PerMallows package is slow because it has to generate the distance frequencies.

# exp_d_ulam_old=function(alpha,n_items){
#   out=PerMallows::expectation.mm(theta=alpha,perm.length=n_items,dist.name="ulam")
#   return(out)
# }

exp_d_ulam <- function(alpha,n_items){ # for n_items<=95
  idx <- 0:(n_items-1)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric=="ulam"][[n_items]]
  norm_const <- exp(get_partition_function(alpha=alpha*n_items,n_items=n_items,
                                        metric="ulam",
                                        cardinalities=card))
  out <- sum(idx*exp(-alpha*idx)*card)/norm_const
  return(out)
}

#/' exp_d_foot
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Footrule distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Footrule metric under the Mallows rank model with the Footrule distance.

exp_d_foot <- function(alpha,n_items){ # for n_items<=50
  idx <- seq(0,floor(n_items^2/2),by=2)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric=="footrule"][[n_items]]
  norm_const <- exp(get_partition_function(alpha=alpha*n_items,n_items=n_items,
                                        metric="footrule",
                                        cardinalities=card))
  # print(length(idx))
  # print(length(card))
  out <- sum(idx*exp(-alpha*idx)*card)/norm_const
  return(out)
}

#/' exp_d_spear
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Spearman distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Spearman metric under the Mallows rank model with the Spearman distance.

exp_d_spear <- function(alpha,n_items){ # for n_items<=14
  idx <- seq(0,2*base::choose(n_items+1,3),by=2)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric=="spearman"][[n_items]]
  norm_const <- exp(get_partition_function(alpha=alpha*n_items,n_items=n_items,
                                        metric="spearman",
                                        cardinalities=card))
  # print(length(idx))
  # print(length(card))
  out <- sum(idx*exp(-alpha*idx)*card)/norm_const
  return(out)
}


# n_items=5
# alpha_vector <- seq(from = 0, to = 10, by = 0.5)
# est=estimate_partition_function(method = "importance_sampling",
#                                 alpha_vector = alpha_vector,
#                                 n_items = n_items, metric = "ulam",
#                                 nmc = 1e3, degree = degree)
#
# n_items=5
# exp(get_partition_function(alpha=10,n_items=5,metric="footrule",cardinalities=ppp$values[[69]]))
# idx=seq(0,floor(n_items^2/2),by=2)
# sum(exp(-10*idx)*ppp$values[[69]])
# sum(exp(-50*idx)*ppp$values[[69]])
# sum(exp(-10/n_items*idx)*ppp$values[[69]])
#
# n_items=5
# exp(get_partition_function(alpha=10,n_items=5,metric="ulam",cardinalities=ppp$values[[5]]))
# idx=0:(n_items-1)
# sum(exp(-10*idx)*ppp$values[[5]])
# sum(exp(-50*idx)*ppp$values[[5]])
# sum(exp(-10/n_items*idx)*ppp$values[[5]])
#
#
# sum(exp(-0.5*idx)*ulam_dist_freq[[n_items]])
# sum(exp(-2.5*idx)*ulam_dist_freq[[n_items]])
#
# temp_n_items=c(1:95,ppp$n_items[51:length(ppp$n_items)])
# temp_metric=c(rep("ulam",95),ppp$metric[51:length(ppp$metric)])
# temp_values=c(ulam_dist_freq[1:95],ppp$values[51:length(ppp$values)])
# temp_type=c(rep("cardinalities",95),ppp$type[51:length(ppp$type)])
# temp_message=c(rep("Using exact partition function.",95),ppp$message[51:length(ppp$message)])
#
# temp=as_tibble(data.frame(n_items=temp_n_items,metric=temp_metric,
#               values=temp_values,type=temp_type,message=temp_message),column_name=names(ppp))
#
# ooo=list(n_items=temp_n_items,metric=temp_metric,
#                values=temp_values,type=temp_type,message=temp_message)
# ppp2=as_tibble(ooo)
