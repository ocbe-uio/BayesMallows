% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_mallows.R
\name{update_mallows}
\alias{update_mallows}
\alias{update_mallows.BayesMallows}
\alias{update_mallows.SMCMallows}
\title{Update a Bayesian Mallows model with new users}
\usage{
update_mallows(model, new_data, ...)

\method{update_mallows}{BayesMallows}(
  model,
  new_data,
  model_options = set_model_options(),
  smc_options = set_smc_options(),
  compute_options = set_compute_options(),
  priors = model$priors,
  ...
)

\method{update_mallows}{SMCMallows}(model, new_data, ...)
}
\arguments{
\item{model}{A model object.}

\item{new_data}{Object returned from \code{\link[=setup_rank_data]{setup_rank_data()}} containing new data.}

\item{...}{Optional arguments. Currently not used.}

\item{model_options}{An object of class \code{"BayesMallowsModelOptions"} returned
from \code{\link[=set_model_options]{set_model_options()}}.}

\item{smc_options}{SMC specific options returned from \code{\link[=set_smc_options]{set_smc_options()}}.}

\item{compute_options}{An object of class \code{"BayesMallowsComputeOptions"}
returned from \code{\link[=set_compute_options]{set_compute_options()}}.}

\item{priors}{An object of class \code{"BayesMallowsPriors"} returned from
\code{\link[=set_priors]{set_priors()}}. Defaults to the priors used in \code{model}.}
}
\value{
An updated model, of class "SMCMallows".
}
\description{
Update a Bayesian Mallows model estimated using the Metropolis-Hastings
algorithm in \code{\link[=compute_mallows]{compute_mallows()}} using the sequential Monte Carlo
algorithm described in
\insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}. This
is useful when new data arrives, and is typically more computationally
efficient than running Metropolis-Hastings from scratch.
}
\seealso{
Other modeling: 
\code{\link{compute_mallows_mixtures}()},
\code{\link{compute_mallows}()}
}
\concept{modeling}