# Specify options for computation

Set parameters related to the Metropolis-Hastings algorithm.

## Usage

``` r
set_compute_options(
  nmc = 2000,
  burnin = NULL,
  alpha_prop_sd = 0.1,
  rho_proposal = c("ls", "swap"),
  leap_size = 1,
  aug_method = c("uniform", "pseudo"),
  pseudo_aug_metric = c("footrule", "spearman"),
  swap_leap = 1,
  alpha_jump = 1,
  aug_thinning = 1,
  clus_thinning = 1,
  rho_thinning = 1,
  include_wcd = FALSE,
  save_aug = FALSE,
  save_ind_clus = FALSE
)
```

## Arguments

- nmc:

  Integer specifying the number of iteration of the Metropolis-Hastings
  algorithm to run. Defaults to `2000`. See
  [`assess_convergence()`](assess_convergence.md) for tools to check
  convergence of the Markov chain.

- burnin:

  Integer defining the number of samples to discard. Defaults to `NULL`,
  which means that burn-in is not set.

- alpha_prop_sd:

  Numeric value specifying the \\\sigma\\ parameter of the lognormal
  proposal distribution used for \\\alpha\\ in the Metropolis-Hastings
  algorithm. The logarithm of the proposed samples will have standard
  deviation given by `alpha_prop_sd`. Defaults to `0.1`.

- rho_proposal:

  Character string specifying the proposal distribution of modal ranking
  \\\rho\\. Defaults to "ls", which means that the leap-and-shift
  algorithm of Vitelli et al. (2018) is used. The other option is
  "swap", which means that the swap proposal of Crispino et al. (2019)
  is used instead.

- leap_size:

  Integer specifying the step size of the distribution defined in
  `rho_proposal` for proposing new latent ranks \\rho\\. Defaults to 1.

- aug_method:

  Augmentation proposal for use with missing data. One of "pseudo" and
  "uniform". Defaults to "uniform", which means that new augmented
  rankings are proposed by sampling uniformly from the set of available
  ranks, see Section 4 in Vitelli et al. (2018) . Setting the argument
  to "pseudo" instead, means that the pseudo-likelihood proposal defined
  in Chapter 5 of Stein (2023) is used instead.

- pseudo_aug_metric:

  String defining the metric to be used in the pseudo-likelihood
  proposal. Only used if `aug_method = "pseudo"`. Can be either
  "footrule" or "spearman", and defaults to "footrule".

- swap_leap:

  Integer specifying the leap size for the swap proposal used for
  proposing latent ranks in the case of non-transitive pairwise
  preference data. Note that leap size for the swap proposal when used
  for proposal the modal ranking \\\rho\\ is given by the `leap_size`
  argument above.

- alpha_jump:

  Integer specifying how many times to sample \\\rho\\ between each
  sampling of \\\alpha\\. In other words, how many times to jump over
  \\\alpha\\ while sampling \\\rho\\, and possibly other parameters like
  augmented ranks \\\tilde{R}\\ or cluster assignments \\z\\. Setting
  `alpha_jump` to a high number can speed up computation time, by
  reducing the number of times the partition function for the Mallows
  model needs to be computed. Defaults to `1`.

- aug_thinning:

  Integer specifying the thinning for saving augmented data. Only used
  when `save_aug = TRUE`. Defaults to `1`.

- clus_thinning:

  Integer specifying the thinning to be applied to cluster assignments
  and cluster probabilities. Defaults to `1`.

- rho_thinning:

  Integer specifying the thinning of `rho` to be performed in the
  Metropolis- Hastings algorithm. Defaults to `1`. `compute_mallows`
  save every `rho_thinning`th value of \\\rho\\.

- include_wcd:

  Logical indicating whether to store the within-cluster distances
  computed during the Metropolis-Hastings algorithm. Defaults to
  `FALSE`. Setting `include_wcd = TRUE` is useful when deciding the
  number of mixture components to include, and is required by
  [`plot_elbow()`](plot_elbow.md).

- save_aug:

  Logical specifying whether or not to save the augmented rankings every
  `aug_thinning`th iteration, for the case of missing data or pairwise
  preferences. Defaults to `FALSE`. Saving augmented data is useful for
  predicting the rankings each assessor would give to the items not yet
  ranked, and is required by [`plot_top_k()`](plot_top_k.md).

- save_ind_clus:

  Whether or not to save the individual cluster probabilities in each
  step. This results in csv files `cluster_probs1.csv`,
  `cluster_probs2.csv`, ..., being saved in the calling directory. This
  option may slow down the code considerably, but is necessary for
  detecting label switching using Stephen's algorithm.

## Value

An object of class `"BayesMallowsComputeOptions"`, to be provided in the
`compute_options` argument to [`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md), or
[`update_mallows()`](update_mallows.md).

## References

Crispino M, Arjas E, Vitelli V, Barrett N, Frigessi A (2019). “A
Bayesian Mallows approach to nontransitive pair comparison data: How
human are sounds?” *The Annals of Applied Statistics*, **13**(1),
492–519. [doi:10.1214/18-aoas1203](https://doi.org/10.1214/18-aoas1203)
.  
  
Stein A (2023). *Sequential Inference with the Mallows Model*. Ph.D.
thesis, Lancaster University.  
  
Vitelli V, Sørensen, Crispino M, Arjas E, Frigessi A (2018).
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research*, **18**(1), 1–49.
<https://jmlr.org/papers/v18/15-481.html>.

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_initial_values()`](set_initial_values.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)
