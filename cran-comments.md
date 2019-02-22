## Resubmission Note
This is a package update. New features include:

* A much improved importance sampling algorithm, which avoids numerical overflow.
* Examples showing how to deal with the label switching problem in mixture models.
* Fixed bug which caused the model with data missing at random to not work properly.
* Fixed bug related to the fact the the Armadillo library uses the C++11 library <random>, which makes the random numbers platform dependent. New implementation uses random number generators from Rcpp, which are platform independent.

## Test Environments
* local OS X install, R 3.5.2
* ubuntu 14.04 on Travis, R 3.5.2
* win-build (devel, release and oldrelease)
* Debian Linux, R-devel, GCC ASAN/UBSAN, via R-hub
* Debian Linux, R-release, GCC, with valgrind, via R-hub

## R CMD CHECK results

* ubuntu 14.04 gave the following NOTE:

  installed size is  8.0Mb
  sub-directories of 1Mb or more:
    libs   7.5Mb

## Downstream Dependencies
There are currently no downstream dependencies for this package.
