## Resubmission Note

A bug in the thinning option for cluster assignment in the case of a single 
cluster has been fixed. I have also reduced the computing time of examples and 
CPU time of tests, as a consequence of pretest errors of a previous submission 
of the same version, see below.

This version has been submitted previously, and got the following pretest errors:

Flavor: r-devel-linux-x86_64-debian-gcc
Check: examples, Result: NOTE
 Examples with CPU (user + system) or elapsed time > 5s
                              user system elapsed
 estimate_partition_function 7.154  0.205   0.932
 Examples with CPU time > 2.5 times elapsed time
                               user system elapsed ratio
 estimate_partition_function  7.154  0.205   0.932 7.896
 compute_mallows_sequentially 4.587  0.100   0.595 7.877
 sample_prior                 2.077  0.059   0.275 7.767

Flavor: r-devel-linux-x86_64-debian-gcc
Check: tests, Result: NOTE
   Running 'testthat.R' [165s/63s]
 Running R code in 'testthat.R' had CPU time 2.6 times elapsed time
 
 


## Test Environments

* windows, devel, release and old-release.
* R-CMD-check via GitHub Actions on windows-latest, macOS-latest, ubuntu-20.04 (release), and ubuntu-20.04 (devel).

## R CMD CHECK results

We got 

0 NOTEs, 0 WARNINGs, 0 ERRORs.


## Downstream Dependencies

The package is reverse imported by 'MSmix' and reverse suggested by 'PlackettLuce'. Running revdep_check() with the 'revdepcheck' package returned OK for both packages.
