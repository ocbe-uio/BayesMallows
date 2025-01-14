# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(BayesMallows)

skip_on_cran_linux <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "false") && .Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Linux") {
    skip("Skipping test on CRAN Linux")
  }
}

test_check("BayesMallows")
