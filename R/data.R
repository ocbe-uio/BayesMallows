#' True ranking of the weights of 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_true_ranking"

#' @title Potato weights assessed visually
#'
#' @description
#' Result of ranking potatoes by weight, where the assessors were only allowed
#' to inspected the potatoes visually. 12 assessors ranked 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_visual"

#' @title Potato weights assessed by hand
#'
#' @description
#' Result of ranking potatoes by weight, where the assessors were
#' allowed to lift the potatoes. 12 assessors ranked 20 potatoes.
#'
#' @family datasets
#' @references \insertRef{liu2019}{BayesMallows}
"potato_weighing"

#' Beach preferences
#'
#' Example dataset from \insertCite{vitelli2018}{BayesMallows}, Section 6.2.
#'
#' @family datasets
#' @references \insertAllCited{}
"beach_preferences"

#' Sushi rankings
#'
#' Complete rankings of 10 types of sushi from 5000 assessors
#' \insertCite{kamishima2003}{BayesMallows}.
#'
#' @family datasets
#' @references \insertAllCited{}
"sushi_rankings"

#' Simulated clustering data
#'
#' Simulated dataset of 60 complete rankings of five items, with three
#' different clusters.
#'
#' @family datasets
"cluster_data"

#' @title Simulated intransitive pairwise preferences
#'
#' @description Simulated dataset based on the [potato_visual] data. Based on
#'   the rankings in [potato_visual], all n-choose-2 = 190 pairs of items were
#'   sampled from each assessor. With probability .9, the pairwise
#'   preference was in agreement with [potato_visual], and with probability .1,
#'   they were in disagreement. Hence, the data generating mechanism was a
#'   Bernoulli error model \insertCite{crispino2019}{BayesMallows} with
#'   \eqn{\theta=0.1}.
#'
#' @family datasets
"bernoulli_data"

#' Sounds data
#'
#' Data from an experiment in which 46 individuals compared 12 different sounds
#' \insertCite{barrett2018}{BayesMallows}. Each assessor was asked multiple
#' times to compare a pair of two sounds, indicating which of the sounds sounded
#' the most like it was human generated. The pairwise preference for each
#' assessor are in general non-transitive. These data inspired the Mallows model
#' for non-transitive pairwise preferences developed by
#' \insertCite{crispino2019}{BayesMallows}.
#'
#' @family datasets
#' @references \insertAllCited{}
"sounds"
