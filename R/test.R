#' Title
#'
#' @param alpha aa
#' @param rho aa
#' @param n_items aa
#' @param partial_ranking aa
#' @param current_ranking aa
#' @param metric aa
#' @param seed aaa
#'
#' @return
#' @export
#'
metropolis_hastings_aug_ranking_wrapper <- function(
  alpha,
  rho,
  n_items,
  partial_ranking,
  current_ranking,
  metric,
  seed = 123L
){
  set.seed(seed)
  metropolis_hastings_aug_ranking(
    alpha,
    rho,
    n_items,
    partial_ranking,
    current_ranking,
    metric
  )
}
