## Resubmission Note
This is a resubmission. On the first submission, we got the following feedback, which is now fixed:

---
  Thanks, we see:
  
     New submission
  
     License components with restrictions and base license permitting such:
       GPL-3 + file LICENSE
     File 'LICENSE':
       GNU General Public License
       ==========================
  
       _Version 3, 29 June 2007_   ..
  
  
  Please only use file LICENSE for additional restrictions for the GPL-3. 
  If there are none, omit it.
  Do not ship the license file itself. It is part of R already.
  
  Best,
  Uwe Ligges
---


## Test Environments
* local OS X install, R 3.5.1
* win-builder (devel and release)

## R CMD CHECK results
The were no ERRORs or WARNINGs in any of the test environments.


On local OS X, there were no NOTEs.

On win-builder (devel and release) there was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Oystein Sorensen <oystein.sorensen.1985@gmail.com>'

New submission

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
