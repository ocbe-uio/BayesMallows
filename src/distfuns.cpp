#include "RcppArmadillo.h"
#include "misc.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Get the distances for computing the partition function given
//' the cardinalities.
//'
//' @param n The number of items.
//' @param cardinalities Number of times each distance appears.
//' @param metric A string. Available options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec get_summation_distances(int n, arma::vec cardinalities,
                                  std::string metric = "footrule") {

  arma::vec distances;

  if(metric == "footrule"){

    if(n > 50) Rcpp::stop("n > 50 currently not supported for footrule");

    int max = floor(pow(n, 2) / 2.0);
    // Sequence from 0 to max with increment 1
    distances = arma::linspace(0, max, max + 1);

  } else if(metric == "spearman") {

    if(n > 13) Rcpp::stop("n > 13 currently not supported for Spearman distance");

    int max = 2 * binomial_coefficient(n, 3);
    distances = arma::linspace(0, max, max + 1);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return distances;
}


//' Compute the distance between two rank vectors.
//'
//' @param r1 A vector of ranks.
//' @param r2 A vector of ranks.
//' @param metric A string. Available options are \code{"footrule"},
//' \code{"kendall"}, \code{"cayley"}, \code{"hamming"} and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar.
//' @details Note that the Spearman distance is the squared L2 norm, whereas
//' the footrule distance is the L1 norm.
//' @keywords internal
//'
// [[Rcpp::export]]
double get_rank_distance(arma::vec r1, arma::vec r2, std::string metric = "footrule"){

  if (r1.n_elem != r2.n_elem){
    Rcpp::stop("r1 and r2 must have the same length");
  }
  int n = r1.n_elem;
  double distance = 0;

  if(metric == "footrule") {

    // Footrule is the sum of absolute values
    distance = arma::norm(r1 - r2, 1);

  } else if(metric == "kendall") {

    // Need loops to compute Kendall distance

    for(int i = 0; i < n; ++i){
      for(int j = 0; j < i; ++j){
        if(((r1(j) > r1(i)) & (r2(j) < r2(i)) ) || ((r1(j) < r1(i)) & (r2(j) > r2(i)))) {
          distance += 1;
        }
      }
    }

  } else if (metric == "cayley") {

    double tmp;

    // This is a C++ translation of Rankcluster::distCayley
    for(int i = 0; i < n; ++i){
      if(r1(i) != r2(i)) {
        distance += 1;
        tmp = r1(i);
        r1(i) = r2(i);
        arma::uvec inds = arma::find(r1 == r2(i));
        r1.elem(inds).fill(tmp);
      }
    }


  } else if (metric == "hamming") {

    return arma::sum(r1 != r2);

  } else if (metric == "spearman") {

    // Spearman distance is the sum of squares
    distance = pow(arma::norm(r1 - r2, 2), 2.0);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return distance;
}

// Compute the distance between all rows in rankings and rho
double rank_dist_matrix(const arma::mat& rankings, const arma::vec& rho, std::string metric){
  int N = rankings.n_cols;
  if(rankings.n_rows != rho.n_elem) Rcpp::stop("rankings and rho have different number of elements");

  double total_distance = 0;

  for(int i = 0; i < N; ++i){
    total_distance += get_rank_distance(rankings.col(i), rho, metric = metric);
  }

  return total_distance;
}


// update the matrix used to speed up clustering
void update_distance_matrix(arma::mat& dist_mat, const arma::mat& rankings, const arma::vec& rho_cluster,
                            const int& n_assessors, const int& cluster_index,
                            const std::string metric){
  for(int i = 0; i < n_assessors; ++i){
    dist_mat(i, cluster_index) = get_rank_distance(rankings.col(i), rho_cluster, metric);
  }

}
