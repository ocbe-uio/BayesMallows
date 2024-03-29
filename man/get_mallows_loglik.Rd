% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mallows_loglik.R
\name{get_mallows_loglik}
\alias{get_mallows_loglik}
\title{Likelihood and log-likelihood evaluation for a Mallows mixture model}
\usage{
get_mallows_loglik(
  rho,
  alpha,
  weights,
  metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
  rankings,
  observation_frequency = NULL,
  log = TRUE
)
}
\arguments{
\item{rho}{A matrix of size \verb{n_clusters x n_items} whose rows are
permutations of the first n_items integers corresponding to the modal
rankings of the Mallows mixture components.}

\item{alpha}{A vector of \code{n_clusters} non-negative scalar specifying the
scale (precision) parameters of the Mallows mixture components.}

\item{weights}{A vector of \code{n_clusters} non-negative scalars specifying
the mixture weights.}

\item{metric}{Character string specifying the distance measure to use.
Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"},
\code{"ulam"}, \code{"footrule"}, and \code{"spearman"}.}

\item{rankings}{A matrix with observed rankings in each row.}

\item{observation_frequency}{A vector of observation frequencies (weights) to apply to
each row in \code{rankings}. This can speed up computation if a large
number of assessors share the same rank pattern. Defaults to \code{NULL},
which means that each row of \code{rankings} is multiplied by 1. If
provided, \code{observation_frequency} must have the same number of elements as there
are rows in \code{rankings}, and \code{rankings} cannot be \code{NULL}.}

\item{log}{A logical; if TRUE, the log-likelihood value is returned,
otherwise its exponential. Default is \code{TRUE}.}
}
\value{
The likelihood or the log-likelihood value corresponding to one or
more observed complete rankings under the Mallows mixture rank model with
distance specified by the \code{metric} argument.
}
\description{
Compute either the likelihood or the log-likelihood value of the
Mallows mixture model parameters for a dataset of complete rankings.
}
\examples{
# Simulate a sample from a Mallows model with the Kendall distance

n_items <- 5
mydata <- sample_mallows(
  n_samples = 100,
  rho0 = 1:n_items,
  alpha0 = 10,
  metric = "kendall")

# Compute the likelihood and log-likelihood values under the true model...
get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = mydata,
  log = FALSE
  )

get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = mydata,
  log = TRUE
  )

# or equivalently, by using the frequency distribution
freq_distr <- compute_observation_frequency(mydata)
get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = freq_distr[, 1:n_items],
  observation_frequency = freq_distr[, n_items + 1],
  log = FALSE
  )

get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = freq_distr[, 1:n_items],
  observation_frequency = freq_distr[, n_items + 1],
  log = TRUE
  )
}
\seealso{
Other rank functions: 
\code{\link{compute_expected_distance}()},
\code{\link{compute_observation_frequency}()},
\code{\link{compute_rank_distance}()},
\code{\link{create_ranking}()},
\code{\link{sample_mallows}()}
}
\concept{rank functions}
