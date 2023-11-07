% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_mallows.R
\name{update_mallows}
\alias{update_mallows}
\alias{update_mallows.BayesMallows}
\alias{update_mallows.SMCMallows}
\title{Update a Bayesian Mallows model with new users}
\usage{
update_mallows(
  model,
  new_rankings,
  n_particles,
  augmentation = "pseudo",
  mcmc_steps = 5
)

\method{update_mallows}{BayesMallows}(
  model,
  new_rankings,
  n_particles,
  augmentation = "pseudo",
  mcmc_steps = 5
)

\method{update_mallows}{SMCMallows}(model, new_rankings)
}
\arguments{
\item{model}{A model object.}

\item{new_rankings}{Matrix containing the new set of observed rankings of size
n_assessors by n_items.}

\item{n_particles}{Integer specifying the number of particles.}

\item{augmentation}{One of "pseudo" and "uniform".}

\item{mcmc_steps}{Number of Metropolis-Hastings steps to apply in sequential
Monte Carlo.}
}
\value{
An updated model, of class "SMCMallows".
}
\description{
Update a Bayesian Mallows model with new users
}
\examples{
# Assume our first dataset contains 100 preferences, observed in 10 batches
data_batch1 <- sushi_rankings[1:100, ]

# Calculating rho and alpha samples
mod1 <- smc_mallows_new_users(
  rankings = data_batch1,
  n_particles = 100,
  timesteps = 10,
  mcmc_steps = 10,
  num_new_obs = 10
)

# We can plot the posterior of alpha
plot(mod1)

# Next, assume a new batch arrives, with 100 new preferences
data_batch2 <- sushi_rankings[101:200, ]

# We update the model in a single timestep. By default, the settings from the
# call to smc_mallows_new_users() is retained.
mod2 <- smc_mallows_update(
  model = mod1,
  rankings = data_batch2
)

plot(mod2)
}
\seealso{
Other modeling: 
\code{\link{compute_mallows_mixtures}()},
\code{\link{compute_mallows}()},
\code{\link{smc_mallows_new_item_rank}()}
}
\concept{modeling}