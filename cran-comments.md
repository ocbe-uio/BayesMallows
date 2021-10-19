## Resubmission Note
This is a resubmission which fixes a critical bug which caused results to be wrong with more than one mixture component in compute_mallows() and compute_mallows_mixtures().

## Test Environments
* local Windows install, R 4.1.1
* windows, win-devel.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

0 NOTE, 0 WARNINGs, 0 ERRORs.

## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
