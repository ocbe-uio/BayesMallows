% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_posterior_intervals.R
\name{compute_posterior_intervals}
\alias{compute_posterior_intervals}
\alias{compute_posterior_intervals.BayesMallows}
\alias{compute_posterior_intervals.SMCMallows}
\title{Compute Posterior Intervals}
\usage{
compute_posterior_intervals(model_fit, ...)

\method{compute_posterior_intervals}{BayesMallows}(
  model_fit,
  parameter = c("alpha", "rho", "cluster_probs"),
  level = 0.95,
  decimals = 3L,
  ...
)

\method{compute_posterior_intervals}{SMCMallows}(
  model_fit,
  parameter = c("alpha", "rho"),
  level = 0.95,
  decimals = 3L,
  ...
)
}
\arguments{
\item{model_fit}{A model object.}

\item{...}{Other arguments. Currently not used.}

\item{parameter}{Character string defining which parameter to compute
posterior intervals for. One of \code{"alpha"}, \code{"rho"}, or
\code{"cluster_probs"}. Default is \code{"alpha"}.}

\item{level}{Decimal number in \eqn{[0,1]} specifying the confidence level.
Defaults to \code{0.95}.}

\item{decimals}{Integer specifying the number of decimals to include in
posterior intervals and the mean and median. Defaults to \code{3}.}
}
\description{
Compute posterior intervals of parameters of interest.
}
\details{
This function computes both the Highest Posterior Density Interval (HPDI),
which may be discontinuous for bimodal distributions, and
the central posterior interval, which is simply defined by the quantiles of the posterior
distribution.
}
\examples{
set.seed(1)
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 3000, burnin = 1000))

# First we compute the interval for alpha
compute_posterior_intervals(model_fit, parameter = "alpha")
# We can reduce the number decimals
compute_posterior_intervals(model_fit, parameter = "alpha", decimals = 2)
# By default, we get a 95 \% interval. We can change that to 99 \%.
compute_posterior_intervals(model_fit, parameter = "alpha", level = 0.99)
# We can also compute the posterior interval for the latent ranks rho
compute_posterior_intervals(model_fit, parameter = "rho")

\dontrun{
  # Posterior intervals of cluster probabilities
  model_fit <- compute_mallows(
    setup_rank_data(sushi_rankings),
    model_options = set_model_options(n_clusters = 5))
  burnin(model_fit) <- 1000

  compute_posterior_intervals(model_fit, parameter = "alpha")

  compute_posterior_intervals(model_fit, parameter = "cluster_probs")
}


}
\references{
\insertAllCited{}
}
\seealso{
Other posterior quantities: 
\code{\link{assign_cluster}()},
\code{\link{compute_consensus}()},
\code{\link{get_acceptance_ratios}()},
\code{\link{heat_plot}()},
\code{\link{plot.BayesMallows}()},
\code{\link{plot.SMCMallows}()},
\code{\link{plot_elbow}()},
\code{\link{plot_top_k}()},
\code{\link{predict_top_k}()},
\code{\link{print.BayesMallows}()}
}
\concept{posterior quantities}
