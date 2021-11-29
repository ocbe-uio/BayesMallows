#' True ranking of the weights of 20 potatoes.
#'
#' @references \insertRef{liu2019}{BayesMallows}
"potato_true_ranking"

#' Result of ranking potatoes by weight, where the assessors were only allowed
#' to inspected the potatoes visually. 12 assessors ranked 20 potatoes.
#'
#' @references \insertRef{liu2019}{BayesMallows}
"potato_visual"

#' Result of ranking potatoes by weight, where the assessors were
#' allowed to lift the potatoes. 12 assessors ranked 20 potatoes.
#'
#' @references \insertRef{liu2019}{BayesMallows}
"potato_weighing"

#' Beach Preferences
#'
#' Example dataset from \insertCite{vitelli2018}{BayesMallows}, Section 6.2.
#'
#' @references \insertAllCited{}
"beach_preferences"

#' Sushi Rankings
#'
#' Complete rankings of 10 types of sushi from
#' 5000 assessors \insertCite{kamishima2003}{BayesMallows}.
#'
#' @references \insertAllCited{}
"sushi_rankings"

#'
#'
#' A synthetic 3D matrix (\code{n_users}, \code{n_items}, \code{Time}) generated
#' using the sample_mallows function. These are test datasets used to run
#' the SMC-Mallows framework for the cases where we know all of the users
#' in our system and their original ranking information are partial rankings.
#' However at some point in time, we observe extra information about
#' an existing user in the form of a rank for an item that was previously
#' not known (\code{NA}). These datasets are very contrived as the first
#' time step (\code{sample_dataset[, , 1]}) we observed the top \code{m / 2}
#' items from each user, where \code{m} is the number of items in a ranking.
#' Then, as we increase the time, we observe the next top ranked item from
#' one user at a time, then the next top ranked item, and so on until we have
#' a complete dataset at \code{sample_dataset[, , Time]}.
#'
#' @references https://github.com/anjastein/SMC-Mallows/tree/main/data
"sample_dataset"
