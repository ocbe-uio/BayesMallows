# Get Acceptance Ratios

Extract acceptance ratio from Metropolis-Hastings algorithm used by
[`compute_mallows()`](compute_mallows.md) or by the move step in
[`update_mallows()`](update_mallows.md) and
[`compute_mallows_sequentially()`](compute_mallows_sequentially.md). If
burnin is not set in the call to
[`compute_mallows()`](compute_mallows.md), the acceptance ratio for all
iterations will be reported. Otherwise the post burnin acceptance ratio
is reported. For the SMC method the acceptance ratios apply to all
iterations, since no burnin is needed in here.

## Usage

``` r
get_acceptance_ratios(model_fit, ...)

# S3 method for class 'BayesMallows'
get_acceptance_ratios(model_fit, ...)

# S3 method for class 'SMCMallows'
get_acceptance_ratios(model_fit, ...)
```

## Arguments

- model_fit:

  A model fit.

- ...:

  Other arguments passed on to other methods. Currently not used.

## Value

A list with elements `alpha_acceptance`, `rho_acceptance`, and
`aug_acceptance`. Each element contains acceptance ratios (between 0
and 1) for the corresponding parameter proposals in the
Metropolis-Hastings algorithm. For models with multiple chains, each
element is a list with one acceptance ratio per chain. Higher values
indicate higher acceptance rates for the Metropolis-Hastings proposals.

## See also

Other posterior quantities: [`assign_cluster()`](assign_cluster.md),
[`compute_consensus()`](compute_consensus.md),
[`compute_posterior_intervals()`](compute_posterior_intervals.md),
[`heat_plot()`](heat_plot.md),
[`plot.BayesMallows()`](plot.BayesMallows.md),
[`plot.SMCMallows()`](plot.SMCMallows.md),
[`plot_elbow()`](plot_elbow.md), [`plot_top_k()`](plot_top_k.md),
[`predict_top_k()`](predict_top_k.md),
[`print.BayesMallows()`](print.BayesMallows.md)

## Examples

``` r
set.seed(1)
mod <- compute_mallows(
  data = setup_rank_data(potato_visual),
  compute_options = set_compute_options(burnin = 200)
)

get_acceptance_ratios(mod)
#> $alpha_acceptance
#> $alpha_acceptance[[1]]
#> [1] 0.7038889
#> 
#> 
#> $rho_acceptance
#> $rho_acceptance[[1]]
#> [1] 0.4716667
#> 
#> 
#> $aug_acceptance
#> $aug_acceptance[[1]]
#> [1] NaN
#> 
#> 
```
