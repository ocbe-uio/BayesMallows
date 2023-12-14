#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Asymptotic Approximation of Partition Function
//'
//' Compute the asymptotic approximation of the logarithm of the partition function,
//' using the iteration algorithm of \insertCite{mukherjee2016;textual}{BayesMallows}.
//'
//' @param alpha_vector A numeric vector of alpha values.
//' @param n_items Integer specifying the number of items.
//' @param metric One of \code{"footrule"} and \code{"spearman"}.
//' @param K Integer.
//' @param n_iterations Integer specifying the number of iterations.
//' @param tol Stopping criterion for algorithm. The previous matrix is subtracted
//' from the updated, and if the maximum absolute relative difference is below \code{tol},
//' the iteration stops.
//'
//' @return A vector, containing the partition function at each value of alpha.
//' @keywords internal
//'
//' @references \insertAllCited{}
//'
// [[Rcpp::export]]
arma::vec asymptotic_partition_function(arma::vec alpha_vector, int n_items, std::string metric,
                                        int K, int n_iterations = 1000, double tol = 1e-9){
  // IPFP procedure
  // Initialize a square matrix where each row/column sums to one
  mat A = ones<mat>(K, K) * 1.0 / K;

  // accu sums all elements of the tensor
  // the % operator computes the element-wise product
  double Z0lim = -2 * std::log(static_cast<double>(K)) - accu(A % log(A));

  // Helper matrix
  mat B(K, K);
  for(int i = 0; i < K; ++i){
    for(int j = 0; j < K; ++j){
      if(metric == "footrule"){
        // Because 0 is the first index, this is really (i + 1) - (j + 1) = i - j
        B(i, j) = -std::abs(static_cast<double>(i - j));
      } else if(metric == "spearman"){
        B(i, j) = -std::pow(static_cast<double>(i - j), 2.0);
      }
    }
  }
  // Divide each element by K
  B = B * 1.0 / K;

  int n_alpha = alpha_vector.n_elem;
  vec result(n_alpha);

  for(int i = 0; i < n_alpha; ++i){
    double alpha = alpha_vector(i);

    A = exp(alpha * B);
    mat A_old = A;
    for(int i = 0; i < n_iterations; ++i){
      // Note: We can use 1-norm because the exponential never gets negative
      // Normalize rows
      A = normalise(A, 1, 1);
      // Normalize columns
      A = normalise(A, 1, 0);

      double diff = abs((A - A_old)/A_old).max();
      if(diff < tol) break;

      A_old = A;
    }

    double Zlim = alpha * accu(B % A) - 2 * std::log(K) - accu(A % log(A));
    double Z0 = sum(log(regspace(1, n_items)));

    result(i) = (Zlim - Z0lim)/K * n_items + Z0;
  }

  return(result);
}
