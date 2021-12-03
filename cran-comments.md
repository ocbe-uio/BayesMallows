## Resubmission Note
This is a resubmission which introduces new sequential Monte Carlo algorithms, describing in a new vignette. In addition, several dependencies have been removed.

## Test Environments
* local Windows install, R 4.1.2
* windows, win-devel.
* Apple Silicon (M1) via rhub.
* valgrind and GCC-UBSAN via rhub.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

0 NOTE, 0 WARNINGs, 0 ERRORs.

## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
