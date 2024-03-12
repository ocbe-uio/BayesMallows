## Resubmission Note

This is a resubmission, fixing a failing unit test on noLD and containing a large number of new features.


## Test Environments

* Local Ubuntu 22.04, R 4.3.2, running R CMD check with --use-valgrind option.
* Local Ubuntu 23.04 with R 4.3.2 built from source with option "--with-valgrind-instrumentation=2", running R CMD check with --use-valgrind option.
* r-devel-san via rocker/r-devel-san, running R CMD check with --use-valgrind option.
* r-devel-san via rocker/r-devel-san.
* local Windows install, R 4.3.2.
* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).
* M1 builder.

## R CMD CHECK results

We got 

0 NOTEs, 0 WARNINGs, 0 ERRORs.


## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
