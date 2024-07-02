#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Dynamic Programming C++ implementation
// of LCS problem

using namespace std;

int longestCommonSubsequence(const arma::vec& r1, const arma::vec& r2)
{
  int n = r1.size();
  int m = r2.size();

  arma::vec prev = arma::zeros(m + 1);
  arma::vec cur = arma::zeros(m + 1);

  for (int idx1 = 1; idx1 < n + 1; idx1++) {
    for (int idx2 = 1; idx2 < m + 1; idx2++) {
      if (r1(idx1 - 1) == r2(idx2 - 1))
        cur(idx2) = 1 + prev(idx2 - 1);
      else
        cur(idx2) = 0 + std::max(cur(idx2 - 1), prev(idx2));
    }
    prev = cur;
  }

  return cur[m];
}

// [[Rcpp::export]]
int test(arma::vec r1, arma::vec r2)
{
  return longestCommonSubsequence(r1, r2);
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
test(c(5,1,4,3,2), c(3,1,2,4,5))
*/
