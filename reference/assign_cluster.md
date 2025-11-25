# Assign Assessors to Clusters

Assign assessors to clusters by finding the cluster with highest
posterior probability.

## Usage

``` r
assign_cluster(model_fit, soft = TRUE, expand = FALSE)
```

## Arguments

- model_fit:

  An object of type `BayesMallows`, returned from
  [`compute_mallows()`](compute_mallows.md).

- soft:

  A logical specifying whether to perform soft or hard clustering. If
  `soft=TRUE`, all cluster probabilities are returned, whereas if
  `soft=FALSE`, only the maximum a posterior (MAP) cluster probability
  is returned, per assessor. In the case of a tie between two or more
  cluster assignments, a random cluster is taken as MAP estimate.

- expand:

  A logical specifying whether or not to expand the rowset of each
  assessor to also include clusters for which the assessor has 0 a
  posterior assignment probability. Only used when `soft = TRUE`.
  Defaults to `FALSE`.

## Value

A dataframe. If `soft = FALSE`, it has one row per assessor, and columns
`assessor`, `probability` and `map_cluster`. If `soft = TRUE`, it has
`n_cluster` rows per assessor, and the additional column `cluster`.

## See also

Other posterior quantities:
[`compute_consensus()`](compute_consensus.md),
[`compute_posterior_intervals()`](compute_posterior_intervals.md),
[`get_acceptance_ratios()`](get_acceptance_ratios.md),
[`heat_plot()`](heat_plot.md),
[`plot.BayesMallows()`](plot.BayesMallows.md),
[`plot.SMCMallows()`](plot.SMCMallows.md),
[`plot_elbow()`](plot_elbow.md), [`plot_top_k()`](plot_top_k.md),
[`predict_top_k()`](predict_top_k.md),
[`print.BayesMallows()`](print.BayesMallows.md)

## Examples

``` r
# Fit a model with three clusters to the simulated example data
set.seed(1)
mixture_model <- compute_mallows(
  data = setup_rank_data(cluster_data),
  model_options = set_model_options(n_clusters = 3),
  compute_options = set_compute_options(nmc = 5000, burnin = 1000)
)

head(assign_cluster(mixture_model))
#>   assessor   cluster probability map_cluster
#> 1        1 Cluster 1     0.20525   Cluster 2
#> 2        1 Cluster 2     0.77775   Cluster 2
#> 3        1 Cluster 3     0.01700   Cluster 2
#> 4        2 Cluster 1     0.10625   Cluster 2
#> 5        2 Cluster 2     0.87125   Cluster 2
#> 6        2 Cluster 3     0.02250   Cluster 2
head(assign_cluster(mixture_model, soft = FALSE))
#>    assessor probability map_cluster
#> 2         1     0.77775   Cluster 2
#> 5         2     0.87125   Cluster 2
#> 8         3     0.91450   Cluster 2
#> 11        4     0.92000   Cluster 2
#> 14        5     0.77675   Cluster 2
#> 17        6     0.91525   Cluster 2
```
