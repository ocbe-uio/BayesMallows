#' True ranking of the weights of 20 potatoes.
#'
#'
#' @format An integer vector with 20 elements
#' @references Insert reference to review paper where the potato data is presented.
"potato_true_ranking"

#' Result of ranking potatoes by weight, where the assessors were only allowed
#' to inspected the potatoes visually. 12 assessors ranked 20 potatoes.
#'
#'
#' @format An integer matrix with 12 rows and 20 columns.
#' @references Insert reference to review paper where the potato data is
#'   presented.
"potato_visual"

#' Result of ranking potatoes by weight, where the assessors were allowed to lift the potatoes. 12 assessors ranked 20 potatoes.
#'
#'
#' @format An integer matrix with 12 rows and 20 columns.
#' @references Insert reference to review paper where the potato data is presented.
"potato_weighing"


#' Importance sampling fits in the form of regression coefficients based on
#' \code{log(Z) = b0 + b1 alpha + ... + b10 alpha^10}. These are to be supplied
#' to \code{\link{run_mcmc}} internally, when the importance sampling estimates are to
#' be used.
"importance_sampling_fit"
