# Get transitive closure

A simple method for showing any transitive closure computed by
[`setup_rank_data()`](setup_rank_data.md).

## Usage

``` r
get_transitive_closure(rank_data)
```

## Arguments

- rank_data:

  An object of class `"BayesMallowsData"` returned from
  [setup_rank_data](setup_rank_data.md).

## Value

A dataframe with transitive closure, if there is any.

## See also

Other preprocessing: [`set_compute_options()`](set_compute_options.md),
[`set_initial_values()`](set_initial_values.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)

## Examples

``` r
# Original beach preferences
head(beach_preferences)
#>   assessor bottom_item top_item
#> 1        1           2       15
#> 2        1           5        3
#> 3        1          13        3
#> 4        1           4        7
#> 5        1           5       15
#> 6        1          12        6
dim(beach_preferences)
#> [1] 1442    3
# We then create a rank data object
dat <- setup_rank_data(preferences = beach_preferences)
# The transitive closure contains additional filled-in preferences implied
# by the stated preferences.
head(get_transitive_closure(dat))
#>   assessor bottom_item top_item
#> 1        1           4        1
#> 2        1           5        1
#> 3        1           7        1
#> 4        1           4        3
#> 5        1           5        3
#> 6        1           7        3
dim(get_transitive_closure(dat))
#> [1] 2921    3
```
