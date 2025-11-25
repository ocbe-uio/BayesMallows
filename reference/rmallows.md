# Sample from the Mallows distribution.

Sample from the Mallows distribution with arbitrary distance metric
using a Metropolis-Hastings algorithm.

## Usage

``` r
rmallows(
  rho0,
  alpha0,
  n_samples,
  burnin,
  thinning,
  leap_size = 1L,
  metric = "footrule"
)
```

## Arguments

- rho0:

  Vector specifying the latent consensus ranking.

- alpha0:

  Scalar specifying the scale parameter.

- n_samples:

  Integer specifying the number of random samples to generate.

- burnin:

  Integer specifying the number of iterations to discard as burn-in.

- thinning:

  Integer specifying the number of MCMC iterations to perform between
  each time a random rank vector is sampled.

- leap_size:

  Integer specifying the step size of the leap-and-shift proposal
  distribution.

- metric:

  Character string specifying the distance measure to use. Available
  options are `"footrule"` (default), `"spearman"`, `"cayley"`,
  `"hamming"`, `"kendall"`, and `"ulam"`.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.
