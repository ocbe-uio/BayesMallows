## Resubmission Note

This is a resubmission containing a large number of new features.


## Test Environments

* Local Ubuntu 23.04 with R 4.3.2 built from source with option "--with-valgrind-instrumentation=2", running R CMD check with --use-valgrind option.
* r-devel-san via rocker/r-devel-san.
* local Windows install, R 4.3.3.
* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).
* M1 builder.

## R CMD CHECK results

We got 

0 NOTEs, 0 WARNINGs, 0 ERRORs.


## Downstream Dependencies

The package is reverse imported by 'MSmix' and reverse suggested by 'PlackettLuce'. Running revdep_check() with the 'revdepcheck' package returned OK for both packages.
