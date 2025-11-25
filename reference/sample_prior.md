# Sample from prior distribution

Function to obtain samples from the prior distributions of the Bayesian
Mallows model. Intended to be given to
[`update_mallows()`](update_mallows.md).

## Usage

``` r
sample_prior(n, n_items, priors = set_priors())
```

## Arguments

- n:

  An integer specifying the number of samples to take.

- n_items:

  An integer specifying the number of items to be ranked.

- priors:

  An object of class "BayesMallowsPriors" returned from
  [`set_priors()`](set_priors.md).

## Value

An object of class "BayesMallowsPriorSample", containing `n` independent
samples of \\\alpha\\ and \\\rho\\.

## See also

Other modeling: [`burnin()`](burnin.md), `burnin<-()`,
[`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md),
[`compute_mallows_sequentially()`](compute_mallows_sequentially.md),
[`update_mallows()`](update_mallows.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# We can use a collection of particles from the prior distribution as
# initial values for the sequential Monte Carlo algorithm.
# Here we start by drawing 1000 particles from the priors, using default
# parameters.
prior_samples <- sample_prior(1000, ncol(sushi_rankings))
# Next, we provide the prior samples to update_mallws(), together
# with the first five rows of the sushi dataset
model1 <- update_mallows(
  model = prior_samples,
  new_data = setup_rank_data(sushi_rankings[1:5, ]))
plot(model1)

# We keep adding more data
model2 <- update_mallows(
  model = model1,
  new_data = setup_rank_data(sushi_rankings[6:10, ]))
plot(model2)

model3 <- update_mallows(
  model = model2,
  new_data = setup_rank_data(sushi_rankings[11:15, ]))
plot(model3)
} # }
```
