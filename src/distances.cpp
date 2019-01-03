#include "RcppArmadillo.h"
#include "subset.h"
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]

double cayley_distance(const arma::vec& r1, const arma::vec& r2){
  double distance = 0;
  int n = r1.n_elem;
  double tmp1;
  arma::vec tmp2 = r1;

  // This is a C++ translation of Rankcluster::distCayley
  for(int i = 0; i < n; ++i){
    if(tmp2(i) != r2(i)) {
      distance += 1;
      tmp1 = tmp2(i);
      tmp2(i) = r2(i);
      arma::uvec inds = arma::find(tmp2 == r2(i));
      tmp2.elem(inds).fill(tmp1);
    }
  }
  return distance;
}

double footrule_distance(const arma::vec& r1, const arma::vec& r2){
  return arma::norm(r1 - r2, 1);
}

double hamming_distance(const arma::vec& r1, const arma::vec& r2){
  return arma::sum(r1 != r2);
}

double kendall_distance(const arma::vec& r1, const arma::vec& r2){
  double distance = 0;
  int n = r1.n_elem;

  for(int i = 0; i < n; ++i){
    for(int j = 0; j < i; ++j){
      if(((r1(j) > r1(i)) & (r2(j) < r2(i)) ) || ((r1(j) < r1(i)) & (r2(j) > r2(i)))) {
        distance += 1;
      }
    }
  }

  return distance;
}

double spearman_distance(const arma::vec& r1, const arma::vec& r2){
  return std::pow(arma::norm(r1 - r2, 2), 2.0);
}

double ulam_distance (const arma::vec& r1, const arma::vec& r2){

  int N = r1.n_elem;

  arma::ivec a = arma::conv_to<arma::ivec>::from(r1);
  arma::ivec b = arma::conv_to<arma::ivec>::from(r2);

  int *p1 = (int*) calloc(N, sizeof (int));
  int *p2 = (int*) calloc(N, sizeof (int));

  int distance;

  for(int i = 0; i < N; ++i){
    p1[i] = static_cast<int>(arma::as_scalar(a(i)) - 1);
    p2[i] = static_cast<int>(arma::as_scalar(b(i)) - 1);
  }

  distance = perm0_distance ( N, p1, p2 );

  free(p1);
  free(p2);
  return static_cast<double>(distance);
}


//' Compute the distance between two rank vectors.
//'
//' @param r1 A vector of ranks.
//' @param r2 A vector of ranks.
//' @param metric A string. Available options are \code{"footrule"},
//' \code{"kendall"}, \code{"cayley"}, \code{"hamming"}, \code{"spearman"}, and \code{"ulam"}.
//' Defaults to \code{"footrule"}.
//' @return A scalar.
//' @details Note that the Spearman distance is the squared L2 norm, whereas
//' the footrule distance is the L1 norm.
//'
//' The Ulam distance uses the SUBSET library developed by John Burkardt, available at http://people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html.
//'
//' The implementation of Cayley distance is based on a \code{C++} translation of \code{Rankcluster::distCayley} \insertCite{Grimonprez2016}{BayesMallows}.
//'
//'
//' @references \insertAllCited{}
//'
//' @keywords internal
//'
//'
// [[Rcpp::export]]
double get_rank_distance(arma::vec r1, arma::vec r2, std::string metric){

  if (r1.n_elem != r2.n_elem){
    Rcpp::stop("r1 and r2 must have the same length");
  }

  if (metric == "cayley") {
    return cayley_distance(r1, r2);
  } else if(metric == "footrule") {
    return footrule_distance(r1, r2);
  } else if (metric == "hamming") {
    return hamming_distance(r1, r2);
  } else if(metric == "kendall") {
    return kendall_distance(r1, r2);
  } else if (metric == "spearman") {
    return spearman_distance(r1, r2);
  } else if (metric == "ulam") {
    return ulam_distance(r1, r2);
  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }
}


// Compute the distance between all rows in rankings and rho, and return the sum
double rank_dist_sum(const arma::mat& rankings, const arma::vec& rho, const std::string& metric){
  return arma::sum(rank_dist_vec(rankings, rho, metric));
}


// Compute the distance between each assessor's ranking and the consensus
arma::vec rank_dist_vec(const arma::mat& rankings, const arma::vec& rho,
                                 const std::string& metric){

  int n = rankings.n_cols;
  arma::vec result = arma::zeros(n);

  for(int i = 0; i < n; ++i){
    result(i) = get_rank_distance(rankings.col(i), rho, metric);
  }
  return(result);
}

