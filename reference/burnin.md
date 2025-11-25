# See the burnin

See the current burnin value of the model.

## Usage

``` r
burnin(model, ...)

# S3 method for class 'BayesMallows'
burnin(model, ...)

# S3 method for class 'BayesMallowsMixtures'
burnin(model, ...)

# S3 method for class 'SMCMallows'
burnin(model, ...)
```

## Arguments

- model:

  A model object.

- ...:

  Optional arguments passed on to other methods. Currently not used.

## Value

An integer specifying the burnin, if it exists. Otherwise `NULL`.

## See also

Other modeling: `burnin<-()`, [`compute_mallows()`](compute_mallows.md),
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
