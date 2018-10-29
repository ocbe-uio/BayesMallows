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
