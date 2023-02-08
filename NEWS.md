# BayesMallows (development versions)

* Replaced deprecated ggplot2::aes_ function with ggplot2::aes.
* Refactoring of SMC functions (#257)

# BayesMallows 1.2.2

* Added `plot.SMCMallows()` method
* Changed default values and argument order on several SMC functions (see PR #269)
* Modifications to internal C++ code to avoid CRAN NOTEs.

# BayesMallows 1.2.1

* PerMallows package has been removed from Imports because it is at risk of
  being removed from CRAN. This means that for Ulam distance with more than
  95 items, the user will have to compute an importance sampling estimate.
* Refactoring of data augmentation function for SMC Mallows.
* Improved documentation of `sample_dataset`

# BayesMallows 1.2.0

* Fixed a bug which caused assess_convergence() to fail with 'parameter =
  "cluster_probs"'.
* Fixed a bug in smc_mallows_new_users_partial() and
  smc_mallows_new_users_partial_alpha_fixed().
* metropolis_hastings_aug_ranking_pseudo() has been deprecated. Please use
  metropolis_hastings_aug_ranking() instead, with pseudo=TRUE.
* smc_mallows_new_users_partial_alpha_fixed(), smc_mallows_new_users_complete(),
  and smc_mallows_new_users_partial() have been deprecated. Please use
  smc_mallows_new_users() instead, and set the type= argument to "complete",
  "partial", or "partial_alpha_fixed".
* smc_mallows_new_item_rank_alpha_fixed() has been deprecated. Please use
  smc_mallows_new_item_rank() instead, with argument alpha_fixed=TRUE.
* Fixed unexpected behavior in leap-and-shift proposal distribution for SMC Mallows,
  causing the function to propose the current rank vector with nonzero probability.
* BayesMallows no longer depends on 'dplyr'.
* Quite extensive internal refactoring of C++ code.
* Function lik_db_mix has been renamed to get_mallows_loglik. lik_db_mix still
  exists as deprecated.
* When no initial rankings are provided, compute_mallows() and compute_mallows_mixtures()
  no use independent initial rho in each cluster. Previously a single initial
  rho was used for all cluster. This should potentially improve convergence, but
  will lead to different results when n_clusters>=2 for a given random number seed.

# BayesMallows 1.1.2

* Fixed an issue with stats::reshape causing an error on R-oldrel.
* Fixed an issue with checking the class of objects, where we now consistently use inherits().
* Internal C++ fixes to comply with CRAN checks.

# BayesMallows 1.1.1

* Fixed C++ errors leading to CRAN issues.

# BayesMallows 1.1.0

* Major update, introducing a whole new class of methods using sequential Monte Carlo. Also reducing the number of dependencies.

# BayesMallows 1.0.4.9001

* This is a major update, with new functions for estimating the Bayesian Mallows
  model using sequential Monte Carlo. The methods are described in the vignette
  titled "SMC-Mallows Tutorial".

# BayesMallows 1.0.4.9000

* Removed a large number of dependencies by converting to base R code. This will
make the package easier to install across a range of systems, and less vulnerable
to changes in other packages.

# BayesMallows 1.0.4

* Incorporates changes since 1.0.3, and also remove PLMIX from Imports.

# BayesMallows 1.0.3.9001

* Fixed bug which caused plot_top_k to fail when plotting clusters.
* Improved the default value of rel_widths argument to plot_top_k.
* Wrote unit tests to check that the bugs don't appear again.

# BayesMallows 1.0.3.9000

* Fixed bug which caused importance sampling to fail when running in parallel.
* Fixed issue with error message when trying to plot error probability when compute_mallows has not been set up to compute error probability.
* Increased number of unit tests.

# BayesMallows 1.0.3

* Fixed critical bug which caused results to be wrong with more than one mixture component in compute_mallows() and compute_mallows_mixtures(). Thanks to Anja Stein for discovering the bug.

# BayesMallows 1.0.2

* Function generate_initial_ranking() now has two additional options for generating random initial rankings. This can help with convergence problems, by allowing the MCMC algorithm to run from a range of different starting points.

# BayesMallows 1.0.1

* Fixes a bug in lik_db_mix and expected_dist, in which the scaling parameter used a different parametrization than the rest of the package. All functions in the package now use consistent parametrization of the Mallows model, as stated in the vignette.

# Bayes Mallows 1.0.0

* Function for computing likelihood added.
* Options for dealing with missing values added, and documentation now states how missing values are dealt with.
* Function rank_freq_distr added, for computing the frequency distribution of ranking patterns.
* Function rank_distance added, for computing the distance between rankings.
* Function expected_dist added for computing expectation of several metrics under the Mallows model.

# BayesMallows 0.5.0

* Function compute_consensus now includes an option for computing consensus of augmented ranks.

# BayesMallows 0.4.4

* Fixed bug in predict_top_k and plot_top_k when using aug_thinning > 1.

# BayesMallows 0.4.3

* Updated README and vignette.

# BayesMallows 0.4.2

* Updating a unit test to make sure BayesMallows is compatible with dplyr version 1.0.0.

# BayesMallows 0.4.1

* Improvement of plotting functions, as noted below.

# BayesMallows 0.4.0.9002
* plot.BayesMallows and plot_elbow no longer print titles automatically.

# BayesMallows 0.4.0.9001
* assess_convergence no longer prints legends for clusters, as the cluster number is essentially arbitrary.

# BayesMallows 0.4.0.9000
* Added CITATION.
* Updated test of random number seed.

# BayesMallows 0.4.0
* Implements all fixes since version 0.3.1 below.
* Fixed typo on y-axis label of elbow plot.
* Fixed an issue which caused the cluster probabilities to differ across platforms, despite using the same seed. https://stackoverflow.com/questions/54822702

# BayesMallows 0.3.1.9005
* Fixed a bug which caused `compute_mallows` not to work (without giving any errors) when `rankings` contained missing values.
* Fixed a bug which caused `compute_mallows` to fail when `preferences` had integer columns.

# BayesMallows 0.3.1.9004
* Changed the name of `save_individual_cluster_probs` to `save_ind_clus`, to save typing.

# BayesMallows 0.3.1.9003
* Added a user prompt asking if the user really wants to save csv files, when `save_individual_cluster_probs = TRUE` in compute_mallows.
* Added `alpha_max`, the truncation of the exponential prior for `alpha`, as a user option in `compute_mallows`.

# BayesMallows 0.3.1.9002
* Added functionality for checking label switching. See `?label_switching` for more info.

# BayesMallows 0.3.1.9001
* The internal function `compute_importance_sampling_estimate` has been updated to avoid numerical overflow. Previously, importance sampling failed at below 200 items. Now it works way above 10,000 items.

# BayesMallows 0.3.1
* This is an update of some parts of the C++ code, to avoid failing the sanitizer checks clang-UBSAN and gcc-UBSAN.

# BayesMallows 0.3.0
* See all bullet points below, since 0.2.0.

# BayesMallows 0.2.0.9006
* `generate_transitive_closure`, `generate_initial_ranking`, and `generate_constraints` now are able to run in parallel.
* Large changes to the underlying code base which should make it more maintainable but not affect the user.

# BayesMallows 0.2.0.9005
* `estimate_partition_function` now has an option to run in parallel, leading to significant speed-up.

# BayesMallows 0.2.0.9004
* Implemented the Bernoulli error model. Set `error_model = "bernoulli"` in `compute_mallows` in order to use it. Examples will come later.

# BayesMallows 0.2.0.9003
* Added parallelization option to `compute_mallows_mixtures` and added `parallel` to **Suggests** field.

# BayesMallows 0.2.0.9002
* Deprecated functions `compute_cp_consensus` and `compute_map_consensus` have been removed. Use `compute_consensus` instead.

# BayesMallows 0.2.0.9001
* Clusters are now `factor` variables sorted according to the cluster number. Hence, in plot legends, "Cluster 10" comes after "Cluster 9", rather than after "Cluster 1" which it used to do until now, because it was a `character`.
* `plot.BayesMallows` no longer contains print statements which forces display of plots. Instead plots are returned from the function. Using `p <- plot(fit)` hence does no longer display a plot, whereas using `plot(fit)` without assigning it to an object, displays a plot. Until now the plot was always shown for `rho` and `alpha`.

# BayesMallows 0.2.0.9000
* `compute_mallows` and `sample_mallows` now support Ulam distance, with argument `metric = "ulam"`.
* Slimmed down the vignette significantly, in order to avoid clang-UBSAN error caused by running the vignette (which was then again caused by `Rcpp`, cf. [this issue](https://github.com/RcppCore/Rcpp/issues/832)). The long vignette is no longer needed in any case, since all the functions are well documented with executable examples.

# BayesMallows 0.2.0
* New release on CRAN, which contains all the updates in 0.1.1, described below.

# BayesMallows 0.1.1.9009
* `Rankcluster` package has been removed from dependencies.

# BayesMallows 0.1.1.9008
* Fixed bug with Cayley distance. For this distance, the computational shortcut on p. 8 of Vitelli et al. (2018), JMLR, does not work. However, it was still used. Now, Cayley distance is always computed with complete rank vectors.
* Fixed bug in the default argument `leap_size` to `compute_mallows`. It used to be `floor(n_items / 5)`, which evaluates to zero when `n_items <= 4`. Updated it to `max(1L, floor(n_items / 5))`.
* Added Hamming distance (`metric = "hamming"`) as an option to `compute_mallows` and `sample_mallows`.

# BayesMallows 0.1.1.9007
* Updated `generate_initial_ranking`, `generate_transitive_closure`, and `sample_mallows` to avoid errors when package `tibble` version 2.0.0 is released. This update is purely internal.

# BayesMallows 0.1.1.9006
* Objects of class `BayesMallows` and `BayesMallowsMixtures` now have default print functions, hence avoiding excessive amounts of informations printed to the console if the user happens to write the name of such an object and press Return.
* `compute_mallows_mixtures` no longer sets `include_wcd = TRUE` by default. The user can choose this argument.
* `compute_mallows` has a new argument `save_clus`, which can be set to `FALSE` for not saving cluster assignments.

# BayesMallows 0.1.1.9005
* `assess_convergence` now automatically plots mixtures.
* `compute_mallows_mixtures` now returns an object of class `BayesMallowsMixtures`.

# BayesMallows 0.1.1.9004
* `assess_convergence` now adds prefix *Assessor* to plots when `parameter = "Rtilde"`.
* `predict_top_k` is now an exported function. Previously it was internal.

# BayesMallows 0.1.1.9003
* `compute_posterior_intervals` now has default `parameter = "alpha"`. Until now, this argument has had no default.
* Argument `type` to `plot.BayesMallows` and `assess_convergence` has been renamed to `parameter`, to be more consistent.

# BayesMallows 0.1.1.9002
* Argument `save_augment_data` to `compute_mallows` has been renamed to `save_aug`.
* `compute_mallows` fills in implied ranks when an assessor has only one missing rank. This avoids unnecessary augmentation in MCMC.
* `generate_ranking` and `generate_ordering` now work with missing ranks.

# BayesMallows 0.1.1.9001
Argument `cluster_assignment_thinning` to `compute_mallows` has been renamed to `clus_thin`.

# BayesMallows 0.1.1.9000
Change the interface for computing consensus ranking. Now, CP and MAP consensus are both computed with the `compute_consensus` function, with argument `type` equal to either `"CP"` or `"MAP"`.
