# Simulated intransitive pairwise preferences

Simulated dataset based on the [potato_visual](potato_visual.md) data.
Based on the rankings in [potato_visual](potato_visual.md), all
n-choose-2 = 190 pairs of items were sampled from each assessor. With
probability .9, the pairwise preference was in agreement with
[potato_visual](potato_visual.md), and with probability .1, they were in
disagreement. Hence, the data generating mechanism was a Bernoulli error
model (Crispino et al. 2019) with \\\theta=0.1\\.

## Usage

``` r
bernoulli_data
```

## Format

An object of class `data.frame` with 2280 rows and 3 columns.

## See also

Other datasets: [`beach_preferences`](beach_preferences.md),
[`cluster_data`](cluster_data.md),
[`potato_true_ranking`](potato_true_ranking.md),
[`potato_visual`](potato_visual.md),
[`potato_weighing`](potato_weighing.md), [`sounds`](sounds.md),
[`sushi_rankings`](sushi_rankings.md)
