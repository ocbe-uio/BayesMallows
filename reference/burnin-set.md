# Set the burnin

Set or update the burnin of a model computed using Metropolis-Hastings.

## Usage

``` r
burnin(model, ...) <- value

# S3 method for class 'BayesMallows'
burnin(model, ...) <- value

# S3 method for class 'BayesMallowsMixtures'
burnin(model, ...) <- value
```

## Arguments

- model:

  An object of class `BayesMallows` returned from
  [`compute_mallows()`](compute_mallows.md) or an object of class
  `BayesMallowsMixtures` returned from
  [`compute_mallows_mixtures()`](compute_mallows_mixtures.md).

- ...:

  Optional arguments passed on to other methods. Currently not used.

- value:

  An integer specifying the burnin. If `model` is of class
  `BayesMallowsMixtures`, a single value will be assumed to be the
  burnin for each model element. Alternatively, `value` can be specified
  as an integer vector of the same length as `model`, and hence a
  separate burnin can be set for each number of mixture components.

## Value

An object of class `BayesMallows` with burnin set.

## See also

Other modeling: [`burnin()`](burnin.md),
[`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md),
[`compute_mallows_sequentially()`](compute_mallows_sequentially.md),
[`sample_prior()`](sample_prior.md),
[`update_mallows()`](update_mallows.md)

## Examples

``` r
set.seed(445)
mod <- compute_mallows(setup_rank_data(potato_visual))
assess_convergence(mod)

burnin(mod)
#> NULL
burnin(mod) <- 1500
burnin(mod)
#> [1] 1500
plot(mod)

#'
models <- compute_mallows_mixtures(
  data = setup_rank_data(cluster_data),
  n_clusters = 1:3)
burnin(models)
#> [[1]]
#> NULL
#> 
#> [[2]]
#> NULL
#> 
#> [[3]]
#> NULL
#> 
burnin(models) <- 100
burnin(models)
#> [[1]]
#> [1] 100
#> 
#> [[2]]
#> [1] 100
#> 
#> [[3]]
#> [1] 100
#> 
burnin(models) <- c(100, 300, 200)
burnin(models)
#> [[1]]
#> [1] 100
#> 
#> [[2]]
#> [1] 300
#> 
#> [[3]]
#> [1] 200
#> 
```
