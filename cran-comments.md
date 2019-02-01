## Resubmission Note
This is a package update, intended to solve a sanitizer issue which occured after last submission. 

The current version of BayesMallows on CRAN was submitted and uploaded on January 30th 2019. Before that submission, the package had issues with clang-UBSAN and gcc-UBSAN related to the vignette building. I removed all executable code from the vignette, and passed R CMD check locally using valgrind. I therefore claimed during the submission process that these issues should now be fixed.

As it turns out, after the submission on January 30th, a new issue occured with clang-UBSAN and gcc-UBSAN. The issue is the following, from https://www.stats.ox.ac.uk/pub/bdr/memtests/gcc-UBSAN/BayesMallows/tests/testthat.Rout:

  > misc.cpp:11:29: runtime error: signed integer overflow: 13 * 479001600 cannot be represented in type 'int'

This issue was not present prior to the submission on January 30th, and is related to the inclusion of new functionality in the package.

First, I offer my apologies for naively claiming that the sanitizer issues were fixed. However, this time I have been able to reproduce the exact issue using rhub::check_with_sanitizers(). In addition, rhub::check_with_sanitizers() returned some warnings which are now also fixed. The version now submitted has hence passed rhub::check_with_sanitizers() with 0 ERRORs and 0 WARNINGs.

The full build log from rhub is available here: https://builder.r-hub.io/status/BayesMallows_0.3.1.tar.gz-389e0a4715b942d5b5bd00134a8c422b

I am really sorry for wasting your time on this, but since the issue reported on CRAN was possible to reproduce on R-hub, I chose to fix it as soon as possible, to avoid further problems.

Below is a summary of all tests that were run on this version.

## Test Environments
* local OS X install, R 3.5.2
* ubuntu 14.04 on Travis, R 3.5.2
* Debian Linux, R-devel, GCC ASAN/UBSAN, via R-hub
* Debian Linux, R-release, GCC, with valgrind, via R-hub

## R CMD CHECK results

* ubuntu 14.04 gave the following NOTE:

  installed size is  8.0Mb
  sub-directories of 1Mb or more:
    libs   7.5Mb

## Downstream Dependencies
There are currently no downstream dependencies for this package.
