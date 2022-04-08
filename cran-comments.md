## Resubmission Note
This is a resubmission. We submitted an update last week which solved C++ issues on CRAN. From the CRAN checks, we now see that these issues have been solved. However, there is still an additional issue on r-oldrel, relating to the stats::reshape() function. We have been able to reproduce this issue using rhub::check_with_roldrel(). This update fixes this issue, which we also have confirmed again with rhub::check_with_roldrel().

We apologize for sending this resubmission so close to the previous, but hope that this is the right thing to do, given that the package currently has errors on r-oldrel, and that we now fix them. We have also fixed an issue with checking the class of objects, where we now consistently use inherits().


## Test Environments
* local Windows install, R 4.1.3
* windows, win-devel.
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
