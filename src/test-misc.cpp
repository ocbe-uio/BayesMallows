/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include "misc.h"

context("Assorted functions in misc.h") {
  test_that("Combinations of aug_method and metric are validated correctly") {
    // Testing wrong combinations
    expect_error(is_pseudo("wrong", "footrule"));
    expect_error(is_pseudo("pseudolikelihood", "ulam"));
    expect_error(is_pseudo("pseudolikelihood", "kendall"));
    expect_error(is_pseudo("pseudolikelihood", "cayley"));
    expect_error(is_pseudo("pseudolikelihood", "hamming"));

    // Testing correct combinations
    expect_true(is_pseudo("pseudolikelihood", "footrule"));
    expect_true(is_pseudo("pseudolikelihood", "spearman"));
    expect_false(is_pseudo("random", "footrule"));
    expect_false(is_pseudo("random", "spearman"));
    expect_false(is_pseudo("random", "ulam"));
    expect_false(is_pseudo("random", "kendall"));
    expect_false(is_pseudo("random", "cayley"));
    expect_false(is_pseudo("random", "hamming"));
  }

}
