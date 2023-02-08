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
#include "leapandshift.h"
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

context("Test leap and shift C++") {
  uvec indices{};
  double prob_backward{}, prob_forward{};
  vec rho{}, assumption{}, rho_proposal{};
  int leap_size{}, n_items{};

  test_that("leap and shift works") {
    n_items = 5;
    rho = regspace(1, n_items);
    leap_size = 1;
    leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                   rho, leap_size, false);
    assumption = {2, 1, 3, 4, 5};

    expect_true(approx_equal(rho_proposal, assumption, "absdiff", 1e-5));
    expect_true(prob_backward == prob_forward); // for leap_size=1
    expect_true(std::abs(prob_backward - 0.3) < 1e-5);
  }

  test_that("leap and shift works") {
    n_items = 10;
    rho = regspace(1, n_items);
    rho_proposal = rho;
    leap_size = 3;
    leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                   rho, leap_size, false);
    assumption = {1, 2, 3, 4, 5, 9, 6, 7, 8, 10};
    expect_true(approx_equal(rho_proposal, assumption, "absdiff", 1e-5));
    expect_true(std::abs(prob_backward - 0.025) < 1e-5);
    expect_true(std::abs(prob_forward - 0.0166667) < 1e-5);
  }
}
