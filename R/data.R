#' True ranking of the weights of 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_true_ranking"

#' Result of ranking potatoes by weight, where the assessors were only allowed
#' to inspected the potatoes visually. 12 assessors ranked 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_visual"

#' Result of ranking potatoes by weight, where the assessors were
#' allowed to lift the potatoes. 12 assessors ranked 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_weighing"

#' Beach Preferences
#'
#' Example dataset from \insertCite{vitelli2018}{BayesMallows}, Section 6.2.
#'
#' @family datasets
#' @references \insertAllCited{}
"beach_preferences"

#' Sushi Rankings
#'
#' Complete rankings of 10 types of sushi from
#' 5000 assessors \insertCite{kamishima2003}{BayesMallows}.
#'
#' @family datasets
#' @references \insertAllCited{}
"sushi_rankings"

#' Updated partial rankings of potatoes
#'
#' Simulated dataset based on \code{\link{potato_visual}}, in which 50 % of the
#' potatoes are initially unranked. Through ten subsequent timesteps, more
#' ranks are added each time, given complete rankings at time 11.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
#'
"potato_partial"
