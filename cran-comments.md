## Resubmission Note
This is a resubmission which introduces new sequential Monte Carlo algorithms, describing in a new vignette. In addition, several dependencies have been removed.

## Test Environments
* local Windows install, R 4.1.2
* windows, win-devel.
* Apple Silicon (M1) via rhub.
* valgrind and GCC-UBSAN via rhub.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

win-devel returned the following NOTE. We have carefully checked that all these DOIs are valid, and can be reached with "curl -I -L" as described in the documentation for "CRAN URL checks". We hence believe these to be false positives.:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Oystein Sorensen <oystein.sorensen.1985@gmail.com>'

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1214/15-aos1389
    From: man/BayesMallows.Rd
          man/asymptotic_partition_function.Rd
          man/estimate_partition_function.Rd
    Status: 500
    Message: Internal Server Error
  URL: https://doi.org/10.1214/18-aoas1203
    From: man/BayesMallows.Rd
          man/compute_mallows.Rd
          man/generate_constraints.Rd
    Status: 500
    Message: Internal Server Error

Found the following (possibly) invalid DOIs:
  DOI: 10.1214/15-AOS1389
    From: DESCRIPTION
    Status: Internal Server Error
    Message: 500
  DOI: 10.1214/18-AOAS1203
    From: DESCRIPTION
    Status: Internal Server Error
    Message: 500
    

The other platforms returned:    

0 NOTE, 0 WARNINGs, 0 ERRORs.

## Downstream Dependencies
The package has one downstream dependency, the package PlackettLuce. This package is not affected by this change, as it does not use any of its functions.
