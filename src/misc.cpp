#include <Rcpp.h>
using namespace Rcpp;

// Function to compute the factorial
// taken from http://www.cplusplus.com/forum/unices/33379/
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
