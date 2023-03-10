## Resubmission Note
This is a resubmission. It adds a heat plot function, removes use of some soft-deprecated ggplot functions, improves documentation, and does internal refactoring.


## Test Environments
* local Windows install, R 4.2.2.
* windows, devel, release and old-release.
* Apple Silicon (M1) via rhub.
* valgrind and GCC-UBSAN via rhub.
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
