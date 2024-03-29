% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_mallows_mixtures.R
\name{compute_mallows_mixtures}
\alias{compute_mallows_mixtures}
\title{Compute Mixtures of Mallows Models}
\usage{
compute_mallows_mixtures(
  n_clusters,
  data,
  model_options = set_model_options(),
  compute_options = set_compute_options(),
  priors = set_priors(),
  initial_values = set_initial_values(),
  pfun_estimate = NULL,
  progress_report = set_progress_report(),
  cl = NULL
)
}
\arguments{
\item{n_clusters}{Integer vector specifying the number of clusters to use.}

\item{data}{An object of class "BayesMallowsData" returned from
\code{\link[=setup_rank_data]{setup_rank_data()}}.}

\item{model_options}{An object of class "BayesMallowsModelOptions" returned
from \code{\link[=set_model_options]{set_model_options()}}.}

\item{compute_options}{An object of class "BayesMallowsComputeOptions"
returned from \code{\link[=set_compute_options]{set_compute_options()}}.}

\item{priors}{An object of class "BayesMallowsPriors" returned from
\code{\link[=set_priors]{set_priors()}}.}

\item{initial_values}{An object of class "BayesMallowsInitialValues" returned
from \code{\link[=set_initial_values]{set_initial_values()}}.}

\item{pfun_estimate}{Object returned from \code{\link[=estimate_partition_function]{estimate_partition_function()}}.
Defaults to \code{NULL}, and will only be used for footrule, Spearman, or
Ulam distances when the cardinalities are not available, cf.
\code{\link[=get_cardinalities]{get_cardinalities()}}.}

\item{progress_report}{An object of class "BayesMallowsProgressReported"
returned from \code{\link[=set_progress_report]{set_progress_report()}}.}

\item{cl}{Optional cluster returned from \code{\link[parallel:makeCluster]{parallel::makeCluster()}}. If
provided, chains will be run in parallel, one on each node of \code{cl}.}
}
\value{
A list of Mallows models of class \code{BayesMallowsMixtures}, with
one element for each number of mixtures that was computed. This object can
be studied with \code{\link[=plot_elbow]{plot_elbow()}}.
}
\description{
Convenience function for computing Mallows models with varying numbers of
mixtures. This is useful for deciding the number of mixtures to use in the
final model.
}
\details{
The \code{n_clusters} argument to \code{\link[=set_model_options]{set_model_options()}} is ignored
when calling \code{compute_mallows_mixtures}.
}
\examples{
# SIMULATED CLUSTER DATA
set.seed(1)
n_clusters <- seq(from = 1, to = 5)
models <- compute_mallows_mixtures(
  n_clusters = n_clusters, data = setup_rank_data(cluster_data),
  compute_options = set_compute_options(nmc = 2000, include_wcd = TRUE))

# There is good convergence for 1, 2, and 3 cluster, but not for 5.
# Also note that there seems to be label switching around the 7000th iteration
# for the 2-cluster solution.
assess_convergence(models)
# We can create an elbow plot, suggesting that there are three clusters, exactly
# as simulated.
burnin(models) <- 1000
plot_elbow(models)

# We now fit a model with three clusters
mixture_model <- compute_mallows(
  data = setup_rank_data(cluster_data),
  model_options = set_model_options(n_clusters = 3),
  compute_options = set_compute_options(nmc = 2000))

# The trace plot for this model looks good. It seems to converge quickly.
assess_convergence(mixture_model)
# We set the burnin to 500
burnin(mixture_model) <- 500

# We can now look at posterior quantities
# Posterior of scale parameter alpha
plot(mixture_model)
plot(mixture_model, parameter = "rho", items = 4:5)
# There is around 33 \% probability of being in each cluster, in agreemeent
# with the data simulating mechanism
plot(mixture_model, parameter = "cluster_probs")
# We can also look at a cluster assignment plot
plot(mixture_model, parameter = "cluster_assignment")

# DETERMINING THE NUMBER OF CLUSTERS IN THE SUSHI EXAMPLE DATA
\dontrun{
  # Let us look at any number of clusters from 1 to 10
  # We use the convenience function compute_mallows_mixtures
  n_clusters <- seq(from = 1, to = 10)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(sushi_rankings),
    compute_options = set_compute_options(include_wcd = TRUE))
  # models is a list in which each element is an object of class BayesMallows,
  # returned from compute_mallows
  # We can create an elbow plot
  burnin(models) <- 1000
  plot_elbow(models)
  # We then select the number of cluster at a point where this plot has
  # an "elbow", e.g., n_clusters = 5.

  # Having chosen the number of clusters, we can now study the final model
  # Rerun with 5 clusters
  mixture_model <- compute_mallows(
    rankings = sushi_rankings,
    model_options = set_model_options(n_clusters = 5),
    compute_options = set_compute_options(include_wcd = TRUE))
  # Delete the models object to free some memory
  rm(models)
  # Set the burnin
  burnin(mixture_model) <- 1000
  # Plot the posterior distributions of alpha per cluster
  plot(mixture_model)
  # Compute the posterior interval of alpha per cluster
  compute_posterior_intervals(mixture_model, parameter = "alpha")
  # Plot the posterior distributions of cluster probabilities
  plot(mixture_model, parameter = "cluster_probs")
  # Plot the posterior probability of cluster assignment
  plot(mixture_model, parameter = "cluster_assignment")
  # Plot the posterior distribution of "tuna roll" in each cluster
  plot(mixture_model, parameter = "rho", items = "tuna roll")
  # Compute the cluster-wise CP consensus, and show one column per cluster
  cp <- compute_consensus(mixture_model, type = "CP")
  cp$cumprob <- NULL
  stats::reshape(cp, direction = "wide", idvar = "ranking",
                 timevar = "cluster", varying = list(as.character(unique(cp$cluster))))

  # Compute the MAP consensus, and show one column per cluster
  map <- compute_consensus(mixture_model, type = "MAP")
  map$probability <- NULL
  stats::reshape(map, direction = "wide", idvar = "map_ranking",
                 timevar = "cluster", varying = list(as.character(unique(map$cluster))))

  # RUNNING IN PARALLEL
  # Computing Mallows models with different number of mixtures in parallel leads to
  # considerably speedup
  library(parallel)
  cl <- makeCluster(detectCores() - 1)
  n_clusters <- seq(from = 1, to = 10)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters,
    rankings = sushi_rankings,
    compute_options = set_compute_options(include_wcd = TRUE),
    cl = cl)
  stopCluster(cl)
}



}
\seealso{
Other modeling: 
\code{\link{burnin}()},
\code{\link{burnin<-}()},
\code{\link{compute_mallows}()},
\code{\link{compute_mallows_sequentially}()},
\code{\link{sample_prior}()},
\code{\link{update_mallows}()}
}
\concept{modeling}
