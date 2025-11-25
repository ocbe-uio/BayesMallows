# BayesMallows: Bayesian Preference Learning with the Mallows Rank Model

An implementation of the Bayesian version of the Mallows rank model
(Vitelli et al., Journal of Machine Learning Research, 2018
<https://jmlr.org/papers/v18/15-481.html>; Crispino et al., Annals of
Applied Statistics, 2019
[doi:10.1214/18-AOAS1203](https://doi.org/10.1214/18-AOAS1203) ;
Sorensen et al., R Journal, 2020
[doi:10.32614/RJ-2020-026](https://doi.org/10.32614/RJ-2020-026) ;
Stein, PhD Thesis, 2023 <https://eprints.lancs.ac.uk/id/eprint/195759>).
Both Metropolis-Hastings and sequential Monte Carlo algorithms for
estimating the models are available. Cayley, footrule, Hamming, Kendall,
Spearman, and Ulam distances are supported in the models. The rank data
to be analyzed can be in the form of complete rankings, top-k rankings,
partially missing rankings, as well as consistent and inconsistent
pairwise preferences. Several functions for plotting and studying the
posterior distributions of parameters are provided. The package also
provides functions for estimating the partition function (normalizing
constant) of the Mallows rank model, both with the importance sampling
algorithm of Vitelli et al. and asymptotic approximation with the IPFP
algorithm (Mukherjee, Annals of Statistics, 2016
[doi:10.1214/15-AOS1389](https://doi.org/10.1214/15-AOS1389) ).

## References

Sørensen Ø, Crispino M, Liu Q, Vitelli V (2020). “BayesMallows: An R
Package for the Bayesian Mallows Model.” *The R Journal*, **12**(1),
324–342.
[doi:10.32614/RJ-2020-026](https://doi.org/10.32614/RJ-2020-026) .

## See also

Useful links:

- <https://github.com/ocbe-uio/BayesMallows>

- <https://ocbe-uio.github.io/BayesMallows/>

- Report bugs at <https://github.com/ocbe-uio/BayesMallows/issues>

## Author

**Maintainer**: Oystein Sorensen <oystein.sorensen.1985@gmail.com>
([ORCID](https://orcid.org/0000-0003-0724-3542))

Authors:

- Waldir Leoncio <w.l.netto@medisin.uio.no>

- Valeria Vitelli <valeria.vitelli@medisin.uio.no>
  ([ORCID](https://orcid.org/0000-0002-6746-0453))

- Marta Crispino <crispino.marta8@gmail.com>

- Qinghua Liu <qinghual@math.uio.no>

- Cristina Mollica <cristina.mollica@uniroma1.it>

- Luca Tardella

- Anja Stein
