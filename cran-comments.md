## Resubmission Note

This is a resubmission, fixing issues with "clang-UBSAN", "gcc-UBSAN", and "valgrind". We have reproduced the errors using the Docker images rocker/r-devel-san, and made the necessary changes. 

We also fix failing test expectations due to platform dependence of C++ code.


## Test Environments
* local Windows install, R 4.3.2.
* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

We got 

0 NOTEs, 0 WARNINGs, 0 ERRORs.


## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
