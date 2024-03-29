% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/burnin.R
\name{burnin<-}
\alias{burnin<-}
\alias{burnin<-.BayesMallows}
\alias{burnin<-.BayesMallowsMixtures}
\title{Set the burnin}
\usage{
burnin(model, ...) <- value

\method{burnin}{BayesMallows}(model, ...) <- value

\method{burnin}{BayesMallowsMixtures}(model, ...) <- value
}
\arguments{
\item{model}{An object of class \code{BayesMallows} returned from
\code{\link[=compute_mallows]{compute_mallows()}} or an object of class \code{BayesMallowsMixtures} returned
from \code{\link[=compute_mallows_mixtures]{compute_mallows_mixtures()}}.}

\item{...}{Optional arguments passed on to other methods. Currently not used.}

\item{value}{An integer specifying the burnin. If \code{model} is of class
\code{BayesMallowsMixtures}, a single value will be assumed to be the burnin
for each model element. Alternatively, \code{value} can be specified as an
integer vector of the same length as \code{model}, and hence a separate burnin
can be set for each number of mixture components.}
}
\value{
An object of class \code{BayesMallows} with burnin set.
}
\description{
Set or update the burnin of a model
computed using Metropolis-Hastings.
}
\examples{
set.seed(445)
mod <- compute_mallows(setup_rank_data(potato_visual))
assess_convergence(mod)
burnin(mod)
burnin(mod) <- 1500
burnin(mod)
plot(mod)
#'
models <- compute_mallows_mixtures(
  data = setup_rank_data(cluster_data),
  n_clusters = 1:3)
burnin(models)
burnin(models) <- 100
burnin(models)
burnin(models) <- c(100, 300, 200)
burnin(models)
}
\seealso{
Other modeling: 
\code{\link{burnin}()},
\code{\link{compute_mallows}()},
\code{\link{compute_mallows_mixtures}()},
\code{\link{compute_mallows_sequentially}()},
\code{\link{sample_prior}()},
\code{\link{update_mallows}()}
}
\concept{modeling}
