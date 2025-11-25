# Setup rank data

Prepare rank or preference data for further analyses.

## Usage

``` r
setup_rank_data(
  rankings = NULL,
  preferences = NULL,
  user_ids = numeric(),
  observation_frequency = NULL,
  validate_rankings = TRUE,
  na_action = c("augment", "fail", "omit"),
  cl = NULL,
  max_topological_sorts = 1,
  timepoint = NULL,
  n_items = NULL
)
```

## Arguments

- rankings:

  A matrix of ranked items, of size `n_assessors x n_items`. See
  [`create_ranking()`](create_ranking.md) if you have an ordered set of
  items that need to be converted to rankings. If `preferences` is
  provided, `rankings` is an optional initial value of the rankings. If
  `rankings` has column names, these are assumed to be the names of the
  items. `NA` values in rankings are treated as missing data and
  automatically augmented; to change this behavior, see the `na_action`
  argument to [`set_model_options()`](set_model_options.md). A vector
  length `n_items` is silently converted to a matrix of length
  `1 x n_items`, and names (if any), are used as column names.

- preferences:

  A data frame with one row per pairwise comparison, and columns
  `assessor`, `top_item`, and `bottom_item`. Each column contains the
  following:

  - `assessor` is a numeric vector containing the assessor index.

  - `bottom_item` is a numeric vector containing the index of the item
    that was disfavored in each pairwise comparison.

  - `top_item` is a numeric vector containing the index of the item that
    was preferred in each pairwise comparison.

  So if we have two assessors and five items, and assessor 1 prefers
  item 1 to item 2 and item 1 to item 5, while assessor 2 prefers item 3
  to item 5, we have the following `df`:

  |              |                 |              |
  |--------------|-----------------|--------------|
  | **assessor** | **bottom_item** | **top_item** |
  | 1            | 2               | 1            |
  | 1            | 5               | 1            |
  | 2            | 5               | 3            |

- user_ids:

  Optional `numeric` vector of user IDs. Only only used by
  [`update_mallows()`](update_mallows.md). If provided, new data can
  consist of updated partial rankings from users already in the dataset,
  as described in Section 6 of Stein (2023) .

- observation_frequency:

  A vector of observation frequencies (weights) to apply do each row in
  `rankings`. This can speed up computation if a large number of
  assessors share the same rank pattern. Defaults to `NULL`, which means
  that each row of `rankings` is multiplied by 1. If provided,
  `observation_frequency` must have the same number of elements as there
  are rows in `rankings`, and `rankings` cannot be `NULL`. See
  [`compute_observation_frequency()`](compute_observation_frequency.md)
  for a convenience function for computing it.

- validate_rankings:

  Logical specifying whether the rankings provided (or generated from
  `preferences`) should be validated. Defaults to `TRUE`. Turning off
  this check will reduce computing time with a large number of items or
  assessors.

- na_action:

  Character specifying how to deal with `NA` values in the `rankings`
  matrix, if provided. Defaults to `"augment"`, which means that missing
  values are automatically filled in using the Bayesian data
  augmentation scheme described in Vitelli et al. (2018) . The other
  options for this argument are `"fail"`, which means that an error
  message is printed and the algorithm stops if there are `NA`s in
  `rankings`, and `"omit"` which simply deletes rows with `NA`s in them.

- cl:

  Optional computing cluster used for parallelization when generating
  transitive closure based on preferences, returned from
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).
  Defaults to `NULL`.

- max_topological_sorts:

  When preference data are provided, multiple rankings will be
  consistent with the preferences stated by each users. These rankings
  are the topological sorts of the directed acyclic graph corresponding
  to the transitive closure of the preferences. This number defaults to
  one, which means that the algorithm stops when it finds a single
  initial ranking which is compatible with the rankings stated by the
  user. By increasing this number, multiple rankings compatible with the
  pairwise preferences will be generated, and one initial value will be
  sampled from this set.

- timepoint:

  Integer vector specifying the timepoint. Defaults to `NULL`, which
  means that a vector of ones, one for each observation, is generated.
  Used by [`update_mallows()`](update_mallows.md) to identify data with
  a given iteration of the sequential Monte Carlo algorithm. If not
  `NULL`, must contain one integer for each row in `rankings`.

- n_items:

  Integer specifying the number of items. Defaults to `NULL`, which
  means that the number of items is inferred from `rankings` or from
  `preferences`. Setting `n_items` manually can be useful with pairwise
  preference data in the SMC algorithm, i.e., when `rankings` is `NULL`
  and `preferences` is non-`NULL`, and contains a small number of
  pairwise preferences for a subset of users and items.

## Value

An object of class `"BayesMallowsData"`, to be provided in the `data`
argument to [`compute_mallows()`](compute_mallows.md).

## Note

Setting `max_topological_sorts` larger than 1 means that many possible
orderings of each assessor's preferences are generated, and one of them
is picked at random. This can be useful when experiencing convergence
issues, e.g., if the MCMC algorithm does not mix properly.

It is assumed that the items are labeled starting from 1. For example,
if a single comparison of the following form is provided, it is assumed
that there is a total of 30 items (`n_items=30`), and the initial
ranking is a permutation of these 30 items consistent with the
preference 29\<30.

|              |                 |              |
|--------------|-----------------|--------------|
| **assessor** | **bottom_item** | **top_item** |
| 1            | 29              | 30           |

If in reality there are only two items, they should be relabeled to 1
and 2, as follows:

|              |                 |              |
|--------------|-----------------|--------------|
| **assessor** | **bottom_item** | **top_item** |
| 1            | 1               | 2            |

## References

Stein A (2023). *Sequential Inference with the Mallows Model*. Ph.D.
thesis, Lancaster University.  
  
Vitelli V, Sørensen, Crispino M, Arjas E, Frigessi A (2018).
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research*, **18**(1), 1–49.
<https://jmlr.org/papers/v18/15-481.html>.

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_compute_options()`](set_compute_options.md),
[`set_initial_values()`](set_initial_values.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md)
