#' Checking for Label Switching in the Mallows Mixture Model
#'
#' @description Label switching may sometimes be a problem when running mixture models.
#' The algorithm by Stephens \insertCite{Stephens2000}{BayesMallows}, implemented
#' in the \code{label.switching} package \insertCite{Papastamoulis2016}{BayesMallows}, allows
#' assessment of label switching after MCMC. At the moment, this is the only avaiable option
#' in the \code{BayesMallows} package. The Stephens algorithms requires the individual cluster
#' probabilities of each assessor to be saved in each iteration of the MCMC algorithm. As this
#' potentially requires much memory, the current implementation of \code{\link{compute_mallows}}
#' saves these cluster probabilities to a csv file in each iteration. The example below shows how
#' to perform such a check for label switching in practice.
#'
#' Beware that this functionality is under development. Later releases might let the user
#' determine the directory and filenames of the csv files.
#'
#' @name label_switching
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/label_switching_example.R
NULL
