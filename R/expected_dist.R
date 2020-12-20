#' Expected value of metrics under a Mallows rank model
#'
#' @description Compute the expectation of several metrics under the Mallows rank model.
#' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model.
#' @param n_items Integer specifying the number of items.
#' @param metric Character string specifying the distance measure to use. Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"}, \code{"ulam"} for \code{n_items<=95}, \code{"footrule"} for \code{n_items<=50} and \code{"spearman"} for \code{n_items<=14}.
#'
#' @return A scalar providing the expected value of the \code{metric} under the Mallows rank model with distance specified by the \code{metric} argument.
#' @export
#'
#' @example /inst/examples/expected_dist_example.R

expected_dist=function(alpha,n_items,metric){
  if(alpha<0){
    stop("alpha must be a non-negative value")
  }else{
    if(metric=="kendall"){
      out=exp_d_tau(alpha,n_items)
    }
    if(metric=="cayley"){
      out=exp_d_cay(alpha,n_items)
    }
    if(metric=="hamming"){
      out=exp_d_ham(alpha,n_items)
    }
    if(metric%in%c("ulam","footrule","spearman")){
      pfd=BayesMallows:::partition_function_data
      card=pfd$values[pfd$metric==metric][[n_items]]
      out=exp(log_expected_dist(alpha=alpha*n_items,n_items=n_items,cardinalities=card,metric=metric))
    }
  }
  return(out)
}
