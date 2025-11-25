# Set prior parameters for Bayesian Mallows model

Set values related to the prior distributions for the Bayesian Mallows
model.

## Usage

``` r
set_priors(gamma = 1, lambda = 0.001, psi = 10, kappa = c(1, 3))
```

## Arguments

- gamma:

  Strictly positive numeric value specifying the shape parameter of the
  gamma prior distribution of \\\alpha\\. Defaults to `1`, thus
  recovering the exponential prior distribution used by (Vitelli et
  al. 2018) .

- lambda:

  Strictly positive numeric value specifying the rate parameter of the
  gamma prior distribution of \\\alpha\\. Defaults to `0.001`. When
  `n_cluster > 1`, each mixture component \\\alpha\_{c}\\ has the same
  prior distribution.

- psi:

  Positive integer specifying the concentration parameter \\\psi\\ of
  the Dirichlet prior distribution used for the cluster probabilities
  \\\tau\_{1}, \tau\_{2}, \dots, \tau\_{C}\\, where \\C\\ is the value
  of `n_clusters`. Defaults to `10L`. When `n_clusters = 1`, this
  argument is not used.

- kappa:

  Hyperparameters of the truncated Beta prior used for error probability
  \\\theta\\ in the Bernoulli error model. The prior has the form
  \\\pi(\theta) = \theta^{\kappa\_{1}} (1 - \theta)^{\kappa\_{2}}\\.
  Defaults to `c(1, 3)`, which means that the \\\theta\\ is a priori
  expected to be closer to zero than to 0.5. See (Crispino et al. 2019)
  for details.

## Value

An object of class `"BayesMallowsPriors"`, to be provided in the
`priors` argument to [`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md), or
[`update_mallows()`](update_mallows.md).

## References

Crispino M, Arjas E, Vitelli V, Barrett N, Frigessi A (2019). “A
Bayesian Mallows approach to nontransitive pair comparison data: How
human are sounds?” *The Annals of Applied Statistics*, **13**(1),
492–519. [doi:10.1214/18-aoas1203](https://doi.org/10.1214/18-aoas1203)
.  
  
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
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)
