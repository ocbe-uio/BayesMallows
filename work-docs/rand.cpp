#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerVector cpp_rbernoulli(int n, double p, int seed) {

  std::default_random_engine generator(seed);
  std::bernoulli_distribution distribution(p);
  IntegerVector out(n);
  for (std::size_t i = 0; i != n; ++i) {
    out[i] = distribution(generator);
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
get_seed <- function() {
  sample.int(.Machine$integer.max, 1)
}

set.seed(1)
cpp_rbernoulli(6, 0.7, seed = get_seed())
cpp_rbernoulli(6, 0.7, seed = get_seed())
set.seed(2)
cpp_rbernoulli(6, 0.7, seed = get_seed())
cpp_rbernoulli(6, 0.7, seed = get_seed())
*/
