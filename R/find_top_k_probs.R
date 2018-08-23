#
# find_top_k_prob <- function(model_fit, burnin, k){
#   stopifnot(burnin < model_fit$nmc)
#
#   n_samples <- model_fit$nmc - burnin
#
#   # Probability of being top-k for rho
#   # Extract post burn-in rows with value <= k
#   rho <- dplyr::filter(model_fit$rho, .data$iteration > burnin, .data$value <= k)
#   rho <- dplyr::group_by(rho, .data$cluster, .data$item)
#   # Find the probabilities
#   rho <- dplyr::summarise(rho, rho_prob_top_k = dplyr::n()/n_samples)
#   rho <- dplyr::ungroup(rho)
#
#   # Probability of top-k ranking
#   rankings <- dplyr::filter(model_fit$augmented_data, .data$iteration > burnin, .data$value <= k)
#   rankings <- dplyr::group_by(rankings, .data$assessor, .data$item)
#   rankings <- dplyr::summarise(rankings, assessor_prob_top_k = dplyr::n()/n_samples)
#   rankings <- dplyr::ungroup(rankings)
#
#
#
# }
