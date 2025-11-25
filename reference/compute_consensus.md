# Compute Consensus Ranking

Compute the consensus ranking using either cumulative probability (CP)
or maximum a posteriori (MAP) consensus (Vitelli et al. 2018) . For
mixture models, the consensus is given for each mixture. Consensus of
augmented ranks can also be computed for each assessor, by setting
`parameter = "Rtilde"`.

## Usage

``` r
compute_consensus(model_fit, ...)

# S3 method for class 'BayesMallows'
compute_consensus(
  model_fit,
  type = c("CP", "MAP"),
  parameter = c("rho", "Rtilde"),
  assessors = 1L,
  ...
)

# S3 method for class 'SMCMallows'
compute_consensus(model_fit, type = c("CP", "MAP"), parameter = "rho", ...)
```

## Arguments

- model_fit:

  A model fit.

- ...:

  Other arguments passed on to other methods. Currently not used.

- type:

  Character string specifying which consensus to compute. Either `"CP"`
  or `"MAP"`. Defaults to `"CP"`.

- parameter:

  Character string defining the parameter for which to compute the
  consensus. Defaults to `"rho"`. Available options are `"rho"` and
  `"Rtilde"`, with the latter giving consensus rankings for augmented
  ranks.

- assessors:

  When `parameter = "rho"`, this integer vector is used to define the
  assessors for which to compute the augmented ranking. Defaults to
  `1L`, which yields augmented rankings for assessor 1.

## References

Vitelli V, Sørensen, Crispino M, Arjas E, Frigessi A (2018).
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research*, **18**(1), 1–49.
<https://jmlr.org/papers/v18/15-481.html>.

## See also

Other posterior quantities: [`assign_cluster()`](assign_cluster.md),
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
# The example datasets potato_visual and potato_weighing contain complete
# rankings of 20 items, by 12 assessors. We first analyse these using the
# Mallows model:
model_fit <- compute_mallows(setup_rank_data(potato_visual))

# Se the documentation to compute_mallows for how to assess the convergence of
# the algorithm. Having chosen burin = 1000, we compute posterior intervals
burnin(model_fit) <- 1000
# We then compute the CP consensus.
compute_consensus(model_fit, type = "CP")
#>      cluster ranking item cumprob
#> 1  Cluster 1       1  P12   1.000
#> 2  Cluster 1       2  P13   1.000
#> 3  Cluster 1       3   P9   1.000
#> 4  Cluster 1       4  P10   0.981
#> 5  Cluster 1       5  P17   0.754
#> 6  Cluster 1       6   P7   0.977
#> 7  Cluster 1       7  P14   1.000
#> 8  Cluster 1       8  P16   1.000
#> 9  Cluster 1       9   P1   0.508
#> 10 Cluster 1      10   P5   0.524
#> 11 Cluster 1      11  P11   0.913
#> 12 Cluster 1      12  P19   0.990
#> 13 Cluster 1      13  P20   0.610
#> 14 Cluster 1      14  P18   1.000
#> 15 Cluster 1      15   P6   0.990
#> 16 Cluster 1      16   P4   0.659
#> 17 Cluster 1      17   P2   0.693
#> 18 Cluster 1      18  P15   1.000
#> 19 Cluster 1      19   P3   1.000
#> 20 Cluster 1      20   P8   1.000
# And we compute the MAP consensus
compute_consensus(model_fit, type = "MAP")
#>      cluster map_ranking item probability
#> 1  Cluster 1           1  P12       0.038
#> 2  Cluster 1           2  P13       0.038
#> 3  Cluster 1           3   P9       0.038
#> 4  Cluster 1           4  P10       0.038
#> 5  Cluster 1           5   P7       0.038
#> 6  Cluster 1           6  P17       0.038
#> 7  Cluster 1           7  P14       0.038
#> 8  Cluster 1           8  P16       0.038
#> 9  Cluster 1           9   P1       0.038
#> 10 Cluster 1          10   P5       0.038
#> 11 Cluster 1          11  P11       0.038
#> 12 Cluster 1          12  P19       0.038
#> 13 Cluster 1          13  P20       0.038
#> 14 Cluster 1          14  P18       0.038
#> 15 Cluster 1          15   P6       0.038
#> 16 Cluster 1          16   P4       0.038
#> 17 Cluster 1          17   P2       0.038
#> 18 Cluster 1          18  P15       0.038
#> 19 Cluster 1          19   P3       0.038
#> 20 Cluster 1          20   P8       0.038

if (FALSE) { # \dontrun{
  # CLUSTERWISE CONSENSUS
  # We can run a mixture of Mallows models, using the n_clusters argument
  # We use the sushi example data. See the documentation of compute_mallows for
  # a more elaborate example
  model_fit <- compute_mallows(
    setup_rank_data(sushi_rankings),
    model_options = set_model_options(n_clusters = 5))
  # Keeping the burnin at 1000, we can compute the consensus ranking per cluster
  burnin(model_fit) <- 1000
  cp_consensus_df <- compute_consensus(model_fit, type = "CP")
  # We can now make a table which shows the ranking in each cluster:
  cp_consensus_df$cumprob <- NULL
  stats::reshape(cp_consensus_df, direction = "wide", idvar = "ranking",
                 timevar = "cluster",
                 varying = list(sort(unique(cp_consensus_df$cluster))))
} # }

if (FALSE) { # \dontrun{
  # MAP CONSENSUS FOR PAIRWISE PREFENCE DATA
  # We use the example dataset with beach preferences.
  model_fit <- compute_mallows(setup_rank_data(preferences = beach_preferences))
  # We set burnin = 1000
  burnin(model_fit) <- 1000
  # We now compute the MAP consensus
  map_consensus_df <- compute_consensus(model_fit, type = "MAP")

  # CP CONSENSUS FOR AUGMENTED RANKINGS
  # We use the example dataset with beach preferences.
  model_fit <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(save_aug = TRUE, aug_thinning = 2))
  # We set burnin = 1000
  burnin(model_fit) <- 1000
  # We now compute the CP consensus of augmented ranks for assessors 1 and 3
  cp_consensus_df <- compute_consensus(
    model_fit, type = "CP", parameter = "Rtilde", assessors = c(1L, 3L))
  # We can also compute the MAP consensus for assessor 2
  map_consensus_df <- compute_consensus(
    model_fit, type = "MAP", parameter = "Rtilde", assessors = 2L)

  # Caution!
  # With very sparse data or with too few iterations, there may be ties in the
  # MAP consensus. This is illustrated below for the case of only 5 post-burnin
  # iterations. Two MAP rankings are equally likely in this case (and for this
  # seed).
  model_fit <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(
      nmc = 1005, save_aug = TRUE, aug_thinning = 1))
  burnin(model_fit) <- 1000
  compute_consensus(model_fit, type = "MAP", parameter = "Rtilde",
                    assessors = 2L)
} # }
```
