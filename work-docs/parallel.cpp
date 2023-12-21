#include <RcppParallel.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::export]]
double matrixSqrt(Rcpp::NumericMatrix orig) {

  double ret{};
  for(const auto a : orig) {
    ret += std::sqrt(a);
  }

  return ret;
}

// [[Rcpp::depends(RcppParallel)]]

using namespace RcppParallel;

struct SquareRoot : public Worker
{
  const RMatrix<double> input;
  double ret{};

  SquareRoot(const Rcpp::NumericMatrix input, Rcpp::NumericMatrix output)
    : input(input), output(output) {}

  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {

  }
};

/*** R

matrixSqrt(matrix(1:10, nrow = 2))


*/
