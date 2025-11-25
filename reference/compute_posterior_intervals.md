# Compute Posterior Intervals

Compute posterior intervals of parameters of interest.

## Usage

``` r
compute_posterior_intervals(model_fit, ...)

# S3 method for class 'BayesMallows'
compute_posterior_intervals(
  model_fit,
  parameter = c("alpha", "rho", "cluster_probs"),
  level = 0.95,
  decimals = 3L,
  ...
)

# S3 method for class 'SMCMallows'
compute_posterior_intervals(
  model_fit,
  parameter = c("alpha", "rho"),
  level = 0.95,
  decimals = 3L,
  ...
)
```

## Arguments

- model_fit:

  A model object.

- ...:

  Other arguments. Currently not used.

- parameter:

  Character string defining which parameter to compute posterior
  intervals for. One of `"alpha"`, `"rho"`, or `"cluster_probs"`.
  Default is `"alpha"`.

- level:

  Decimal number in \\\[0,1\]\\ specifying the confidence level.
  Defaults to `0.95`.

- decimals:

  Integer specifying the number of decimals to include in posterior
  intervals and the mean and median. Defaults to `3`.

## Details

This function computes both the Highest Posterior Density Interval
(HPDI), which may be discontinuous for bimodal distributions, and the
central posterior interval, which is simply defined by the quantiles of
the posterior distribution.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

Other posterior quantities: [`assign_cluster()`](assign_cluster.md),
[`compute_consensus()`](compute_consensus.md),
[`get_acceptance_ratios()`](get_acceptance_ratios.md),
[`heat_plot()`](heat_plot.md),
[`plot.BayesMallows()`](plot.BayesMallows.md),
[`plot.SMCMallows()`](plot.SMCMallows.md),
[`plot_elbow()`](plot_elbow.md), [`plot_top_k()`](plot_top_k.md),
[`predict_top_k()`](predict_top_k.md),
[`print.BayesMallows()`](print.BayesMallows.md)

## Examples

``` r
set.seed(1)
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 3000, burnin = 1000))

# First we compute the interval for alpha
compute_posterior_intervals(model_fit, parameter = "alpha")
#>   parameter   mean median           hpdi central_interval
#> 1     alpha 10.778 10.818 [9.165,12.566]   [8.929,12.351]
# We can reduce the number decimals
compute_posterior_intervals(model_fit, parameter = "alpha", decimals = 2)
#>   parameter  mean median         hpdi central_interval
#> 1     alpha 10.78  10.82 [9.16,12.57]     [8.93,12.35]
# By default, we get a 95 % interval. We can change that to 99 %.
compute_posterior_intervals(model_fit, parameter = "alpha", level = 0.99)
#>   parameter   mean median           hpdi central_interval
#> 1     alpha 10.778 10.818 [8.276,12.748]   [8.276,12.748]
# We can also compute the posterior interval for the latent ranks rho
compute_posterior_intervals(model_fit, parameter = "rho")
#>    parameter item mean median      hpdi central_interval
#> 1        rho   P1   10     10    [9,12]           [9,12]
#> 2        rho   P2   17     17   [16,18]          [16,18]
#> 3        rho   P3   19     19      [19]             [19]
#> 4        rho   P4   16     16   [16,18]          [16,18]
#> 5        rho   P5    9      9 [3][9,11]           [3,11]
#> 6        rho   P6   15     15      [15]          [15,16]
#> 7        rho   P7    6      6     [5,7]            [5,7]
#> 8        rho   P8   20     20      [20]             [20]
#> 9        rho   P9    3      3     [3,4]            [3,4]
#> 10       rho  P10    4      4  [4,5][7]            [3,7]
#> 11       rho  P11   11     11    [9,12]           [9,12]
#> 12       rho  P12    1      1       [1]              [1]
#> 13       rho  P13    2      2       [2]              [2]
#> 14       rho  P14    7      7     [6,7]            [6,7]
#> 15       rho  P15   18     18   [17,18]          [17,18]
#> 16       rho  P16    8      8     [8,9]            [8,9]
#> 17       rho  P17    5      5  [4,6][8]            [4,8]
#> 18       rho  P18   14     14   [13,14]          [13,14]
#> 19       rho  P19   12     12   [10,12]          [10,12]
#> 20       rho  P20   13     13   [13,14]          [13,14]

if (FALSE) { # \dontrun{
  # Posterior intervals of cluster probabilities
  model_fit <- compute_mallows(
    setup_rank_data(sushi_rankings),
    model_options = set_model_options(n_clusters = 5))
  burnin(model_fit) <- 1000

  compute_posterior_intervals(model_fit, parameter = "alpha")

  compute_posterior_intervals(model_fit, parameter = "cluster_probs")
} # }

```
