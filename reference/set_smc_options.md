# Set SMC compute options

Sets the SMC compute options to be used in
[`update_mallows.BayesMallows()`](update_mallows.md).

## Usage

``` r
set_smc_options(
  n_particles = 1000,
  mcmc_steps = 5,
  resampler = c("stratified", "systematic", "residual", "multinomial"),
  latent_sampling_lag = NA_integer_,
  max_topological_sorts = 1
)
```

## Arguments

- n_particles:

  Integer specifying the number of particles.

- mcmc_steps:

  Number of MCMC steps to be applied in the resample-move step.

- resampler:

  Character string defining the resampling method to use. One of
  "stratified", "systematic", "residual", and "multinomial". Defaults to
  "stratified". While multinomial resampling was used in Stein (2023) ,
  stratified, systematic, or residual resampling typically give lower
  Monte Carlo error (Douc and Cappe 2005; Hol et al. 2006; Naesseth et
  al. 2019) .

- latent_sampling_lag:

  Parameter specifying the number of timesteps to go back when
  resampling the latent ranks in the move step. See Section 6.2.3 of
  (Kantas et al. 2015) for details. The \\L\\ in their notation
  corresponds to `latent_sampling_lag`. See more under Details. Defaults
  to `NA`, which means that all latent ranks from previous timesteps are
  moved. If set to `0`, no move step is applied to the latent ranks.

- max_topological_sorts:

  User when pairwise preference data are provided, and specifies the
  maximum number of topological sorts of the graph corresponding to the
  transitive closure for each user will be used as initial ranks.
  Defaults to 1, which means that all particles get the same initial
  augmented ranking. If larger than 1, the initial augmented ranking for
  each particle will be sampled from a set of maximum size
  `max_topological_sorts`. If the actual number of topological sorts
  consists of fewer rankings, then this determines the upper limit.

## Value

An object of class "SMCOptions".

## Lag parameter in move step

The parameter `latent_sampling_lag` corresponds to \\L\\ in (Kantas et
al. 2015) . Its use in this package is can be explained in terms of
Algorithm 12 in (Stein 2023) . The relevant line of the algorithm is:

**for** \\j = 1 : M\_{t}\\ **do**  
**M-H step:** update \\\tilde{\mathbf{R}}\_{j}^{(i)}\\ with proposal
\\\tilde{\mathbf{R}}\_{j}' \sim q(\tilde{\mathbf{R}}\_{j}^{(i)} \|
\mathbf{R}\_{j}, \boldsymbol{\rho}\_{t}^{(i)}, \alpha\_{t}^{(i)})\\.  
**end**

Let \\L\\ denote the value of `latent_sampling_lag`. With this
parameter, we modify for algorithm so it becomes

**for** \\j = M\_{t-L+1} : M\_{t}\\ **do**  
**M-H step:** update \\\tilde{\mathbf{R}}\_{j}^{(i)}\\ with proposal
\\\tilde{\mathbf{R}}\_{j}' \sim q(\tilde{\mathbf{R}}\_{j}^{(i)} \|
\mathbf{R}\_{j}, \boldsymbol{\rho}\_{t}^{(i)}, \alpha\_{t}^{(i)})\\.  
**end**

This means that with \\L=0\\ no move step is performed on any latent
ranks, whereas \\L=1\\ means that the move step is only applied to the
parameters entering at the given timestep. The default,
`latent_sampling_lag = NA` means that \\L=t\\ at each timestep, and
hence all latent ranks are part of the move step at each timestep.

## References

Douc R, Cappe O (2005). “Comparison of resampling schemes for particle
filtering.” In *ISPA 2005. Proceedings of the 4th International
Symposium on Image and Signal Processing and Analysis, 2005.*.
[doi:10.1109/ispa.2005.195385](https://doi.org/10.1109/ispa.2005.195385)
, <http://dx.doi.org/10.1109/ISPA.2005.195385>.  
  
Hol JD, Schon TB, Gustafsson F (2006). “On Resampling Algorithms for
Particle Filters.” In *2006 IEEE Nonlinear Statistical Signal Processing
Workshop*.
[doi:10.1109/nsspw.2006.4378824](https://doi.org/10.1109/nsspw.2006.4378824)
, <http://dx.doi.org/10.1109/NSSPW.2006.4378824>.  
  
Kantas N, Doucet A, Singh SS, Maciejowski J, Chopin N (2015). “On
Particle Methods for Parameter Estimation in State-Space Models.”
*Statistical Science*, **30**(3). ISSN 0883-4237,
[doi:10.1214/14-sts511](https://doi.org/10.1214/14-sts511) ,
<http://dx.doi.org/10.1214/14-STS511>.  
  
Naesseth CA, Lindsten F, Schön TB (2019). “Elements of Sequential Monte
Carlo.” *Foundations and Trends® in Machine Learning*, **12**(3),
187–306. ISSN 1935-8245,
[doi:10.1561/2200000074](https://doi.org/10.1561/2200000074) ,
<http://dx.doi.org/10.1561/2200000074>.  
  
Stein A (2023). *Sequential Inference with the Mallows Model*. Ph.D.
thesis, Lancaster University.

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_compute_options()`](set_compute_options.md),
[`set_initial_values()`](set_initial_values.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`setup_rank_data()`](setup_rank_data.md)
