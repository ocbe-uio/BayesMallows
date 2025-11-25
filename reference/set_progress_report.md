# Set progress report options for MCMC algorithm

Specify whether progress should be reported, and how often.

## Usage

``` r
set_progress_report(verbose = FALSE, report_interval = 1000)
```

## Arguments

- verbose:

  Boolean specifying whether to report progress or not. Defaults to
  `FALSE`.

- report_interval:

  Strictly positive number specifying how many iterations of MCMC should
  be run between each progress report. Defaults to `1000`.

## Value

An object of class `"BayesMallowsProgressReport"`, to be provided in the
`progress_report` argument to [`compute_mallows()`](compute_mallows.md)
and [`compute_mallows_mixtures()`](compute_mallows_mixtures.md).

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_compute_options()`](set_compute_options.md),
[`set_initial_values()`](set_initial_values.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)
