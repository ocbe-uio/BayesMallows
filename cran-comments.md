## Resubmission Note

This is a resubmission. It fixes a bug related to an argument that was not 
properly forwarded, it allows the sequential Monte Carlo functions to use 
exact partition functions, removes deprecated functions, adds link to 
website, and improves the unit test outputs.


## Test Environments
* local Windows install, R 4.3.1.
* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

We got 

1 NOTE, 0 WARNINGs, 0 ERRORs.

The note was, slightly dependent on platform

* checking installed package size ... NOTE
  installed size is 23.3Mb
  sub-directories of 1Mb or more:
    doc    1.7Mb
    libs  20.8Mb

## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
