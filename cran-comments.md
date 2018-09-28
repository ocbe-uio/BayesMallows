## Comment to Package Reviewer
This is my second submission. My first submission resulted in a pretest 
NOTE on Debian regarding Makevars.win. The issue should now be fixed.

## Test Environments
* local OS X install, R 3.5.1
* win-builder (devel and release)

## R CMD CHECK results
The were no ERRORs or WARNINGs in any of the test environments.


On local OS X, there were no NOTEs.

On win-builder (devel and release) there was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Oystein Sorensen <oystein.sorensen.1985@gmail.com>'
Possibly mis-spelled words in DESCRIPTION:
  Cayley (7:190)
  IPFP (7:748)
  Mukherjee (7:764)
  Vitelli (7:83, 7:695)
  al (7:94, 7:706)
  et (7:91, 7:703)
  footrule (7:198)

## Downstream Dependencies
There are currently no downstream dependencies for this package.
