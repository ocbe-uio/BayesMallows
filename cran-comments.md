## Resubmission Note

An error in compute_mallows_loglik when the number of clusters is more than one has been corrected.

## Test Environments

* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

We got 

0 NOTEs, 0 WARNINGs, 0 ERRORs.


## Downstream Dependencies

The package is reverse imported by 'MSmix' and reverse suggested by 'PlackettLuce'. Running revdep_check() with the 'revdepcheck' package returned OK for both packages.
