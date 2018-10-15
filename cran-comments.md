## Resubmission Note
This is a package update. It attempts to fix the following check results on CRAN:

* `r-devel-linux-x86_64-debian-clang` and `r-devel-linux-x86_64-fedora-clang` currently issue a WARNING. The cause of the warning, on line 127 in `distfuns.cpp`, should now be fixed.
* `r-patched-solaris-x86` currently fails. The `C++` code has been updated according to $1.6.4 of 'Writing R Extensions', and hopefully will succeed with this update.
* `r-oldrel-windows-ix86+x86_64` and `r-oldrel-osx-x86_64` fail because `stats (>= 3.5.0)` was specified under IMPORTS. This has now been changed to `stats`, not explicitly demanding version 3.5.0 or higher.

## Test Environments
* local OS X install, R 3.5.1
* ubuntu 14.04 (on travis-ci), R 3.5.1
* win-builder (devel and release)

## R CMD CHECK results
The were no ERRORs, WARNINGs or NOTEs in any of the test environments.

## Downstream Dependencies
There are currently no downstream dependencies for this package.
