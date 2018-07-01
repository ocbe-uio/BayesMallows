#' BayesMallows: Bayesian Preference Learning with the Mallows Rank Model.
#'
#' @description The BayesMallows package provides functionality for fully
#'   Bayesian analysis of preference or rank data. It provides three categories
#'   of functions:
#'
#'   \itemize{ \item Functions for sampling from the posterior distributions of
#'   the Mallows rank model. \item Functions for studying the posterior
#'   distributions of the Mallows rank model. \item Utility functions for
#'   working with rank data. }
#'
#'   These three classes of functions are described below. Click on the links to
#'   the respective functions to see examples.
#'
#' @section Sampling from the Posterior Distribution of the Mallows Rank Model:
#'   The function \code{\link{compute_mallows}} is the main function for
#'   obtaining samples from the posterior distribution. It is a wrapper to an
#'   underlying C++ function which runs a Metropolis-Hastings algorithm, and
#'   returns an object of class \code{BayesMallows} which contains all the data
#'   describing the posterior distribution. In order to tune the parameters of
#'   this algorithm, it is recommended to start by using the function
#'   \code{\link{assess_convergence}}.
#'
#' @section Studying the Posterior Distribution of the Mallows Rank Model: The
#'   function \code{\link{plot.BayesMallows}} is an S3 method for plotting the
#'   object returned from \code{\link{compute_mallows}}.
#'   \code{summary.BayesMallows} will soon be implemented, returning a dataframe
#'   describing the posterior distribution.
#'
#' @section Utility Functions for Working with Rank Data: The function
#'   \code{\link{get_rank_distance}} computes the distance between two rank
#'   vectors for five different metrics. \code{\link{get_partition_function}}
#'   computes the logarithm of the partition function for the Mallows rank model
#'   for the same five metrics.
#'
#' @section Example Data: The BayesMallows package comes with example data from
#'   an experiment in which the participants were asked to rank a set of
#'   potatoes based on their weight. The datasets are
#'   \code{\link{potato_weighing}}, \code{\link{potato_visual}}, and the true
#'   weights are in the vector \code{\link{potato_true_ranking}}. INSERT
#'   REFERENCE
#'
#' @references \insertRef{vitelli2018}{BayesMallows}
#'
#' @docType package
#' @name BayesMallows
NULL
