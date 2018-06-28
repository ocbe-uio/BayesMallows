#include "RcppArmadillo.h"
#include "math.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// We put all the function declarations here for clarity
int binomial_coefficient(int n, int k);
double get_rank_distance(arma::rowvec, arma::rowvec, std::string);
double get_partition_function(double, Rcpp::List, std::string);
Rcpp::List run_mcmc(arma::Mat<int>, arma::vec, std::string);
arma::vec get_summation_distances(int, arma::vec, std::string);
Rcpp::List leap_and_shift(arma::rowvec, int);
double rank_dist_matrix(arma::mat, arma::rowvec, std::string);

//' Worker function for computing the posterior distribtuion.
//'
//' @param R A set of complete rankings.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function.
//' TODO: Make this argument optional.
//' @param metric The distance metric to use.
//' @param nmc Number of Monte Carlo samples.
//' @param L Leap-and-shift step size.
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat R, arma::vec cardinalities,
                    std::string metric = "footrule", int nmc = 10,
                    int L = 1){

  // The number of items ranked
  int n = R.n_cols;

  // First we find the vector of distances used to compute the partition function
  arma::vec distances = get_summation_distances(n, cardinalities, metric);

  // Declare the matrix to hold the latent ranks
  arma::mat rho = arma::zeros<arma::mat>(nmc, n);

  // Set the initial latent rank value
  rho.row(0) = arma::linspace<arma::rowvec>(1, n, n);

  // Declare the vector to hold the scaling parameter alpha
  arma::vec alpha = arma::zeros<arma::vec>(nmc);

  // Set the initial alpha value
  alpha(0) = 1;

  // Number of times alpha and rho proposals are accepted
  int acceptance_rho = 0, acceptance_alpha = 0;

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0
  for(int t = 1; t < nmc; ++t){
    // Save current parameter values
    arma::rowvec rho_old = rho.row(t - 1);
    double alpha_old = alpha(t - 1);

    // Sample a rank proposal
    Rcpp::List ls_proposal = leap_and_shift(rho_old, L);

    // Save some of the variables
    arma::rowvec rho_proposal = ls_proposal["proposal"];
    arma::uvec indices = ls_proposal["indices"];
    double prob_backward = ls_proposal["prob_backward"];
    double prob_forward = ls_proposal["prob_forward"];

    // Compute the distances to current and proposed ranks
    double dist_new = rank_dist_matrix(R.cols(indices), rho_proposal(indices), metric);
    double dist_old = rank_dist_matrix(R.cols(indices), rho_old(indices), metric);

    // Metropolis-Hastings ratio
    double ratio = - alpha_old / n * (dist_new - dist_old) + log(prob_backward) - log(prob_forward);

    // Draw a uniform random number
    double u = log(arma::randu<double>());

    if(ratio > u){
      rho.row(t) = rho_proposal;
      acceptance_rho += 1;
    } else {
      rho.row(t) = rho.row(t - 1);
    }

    alpha(t) = alpha(t - 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("alpha") = alpha
  );
}

//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param alpha The value of the alpha parameter.
//' @param summation_sequences List with elements \code{distances} and
//' \code{cardinalities}, both of type arma::vec.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
// [[Rcpp::export]]
double get_partition_function(double alpha, Rcpp::List summation_sequences,
                              std::string metric = "footrule"){



  double log_z_n = 0;

  if(metric == "footrule") {
    log_z_n = 1;


  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return log_z_n;
}

//' Get the distances for computing the partition function given
//' the cardinalities.
//'
//' @param alpha The value of the alpha parameter.
//' @param summation_sequences List with elements \code{distances} and
//' \code{cardinalities}, both of type arma::vec.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
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
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar.
//' @details Note that the Spearman distance is the squared L2 norm, whereas
//' the footrule distance is the L1 norm.
//' @export
// [[Rcpp::export]]
double get_rank_distance(arma::rowvec r1, arma::rowvec r2, std::string metric = "footrule"){

  if (r1.n_elem != r2.n_elem){
    Rcpp::Rcout << "r1 and r2 must have the same length" << std::endl;
    exit(1);
  }
  int n = r1.n_elem;
  double distance = 0;

  if(metric == "footrule") {

    // Footrule is the one-norm
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

  } else if (metric == "spearman") {

    // Spearman distance is the squared L2 norm
    distance = pow(arma::norm(r1 - r2, 2), 2);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return distance;
}




int binomial_coefficient(int n, int k){
  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k )
    k = n - k;

  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}



Rcpp::List leap_and_shift(arma::rowvec rho, int L){

  // Declare the proposed rank vector
  arma::rowvec proposal = rho;

  // Help vectors
  arma::vec support, indices;

  // Number of items
  int n = rho.n_elem;

  // Other helper variables
  int u, index;
  double delta_r, prob_forward, prob_backward, support_new;

  // Leap 1
  // 1, sample u randomly between 1 and n
  u = arma::as_scalar(arma::randi(1, arma::distr_param(1, n)));

  // 2, compute the set S for sampling the new rank
  // Defining versions of L and n converted to double, to avoid duplication in code
  double dobL = static_cast<double>(L);
  double dobn = static_cast<double>(n);

  // Defining linspace lengths here to avoid duplication in code
  double length1 = std::min(rho(u - 1) - 1, dobL);
  double length2 = std::min(n - rho(u - 1), dobL);

  if((rho(u - 1) > 1) & (rho(u - 1) < n)){
    support = arma::join_cols(
      arma::linspace(
        std::max(1.0, rho(u - 1) - L), rho(u - 1) - 1, length1
    ),
    arma::linspace(
      rho(u - 1) + 1, std::min(dobn, rho(u - 1) + L), length2
    )
    );
  } else if(rho(u - 1) == 1){
    support = arma::linspace(
      rho(u - 1) + 1,
      std::min(dobn, rho(u - 1) + L),
      length2
    );

  } else if(rho(u - 1) == n){
    support = arma::linspace(
      std::max(1.0, rho(u - 1) - L), rho(u - 1) - 1,
      length1);
  }

  // 3. assign a random element of the support set, this completes the leap step
  index = arma::as_scalar(arma::randi(1, arma::distr_param(0, support.n_elem-1)));
  // Picked element index-1 from the support set
  proposal(u-1) = support(index);

  // Compute the associated transition probabilities (BEFORE THE SHIFT STEP, WHICH IS DETERMINISTIC --> EASIER)
  if(std::abs(proposal(u - 1) - rho(u - 1)) == 1){
    // in this case the transition probabilities coincide! (and in fact for L = 1 the L&S is symmetric)
    support_new = std::min(proposal(u - 1) - 1, dobL) + std::min(n - proposal(u - 1), dobL);
    prob_forward = 1.0 / (n * support.n_elem) + 1.0 / (n * support_new);
    prob_backward = prob_forward;
  } else {
    // P(proposed|current)
    prob_forward = 1.0 / (n * support.n_elem);
    // P(current|proposed)
    support_new = std::min(proposal(u - 1) - 1, dobL) + std::min(n - proposal(u-1), dobL);
    prob_backward = 1.0 / (n * support_new);
  }

  // Shift step:
  delta_r = proposal(u - 1) - rho(u - 1);
  indices = arma::zeros(std::abs(delta_r) + 1);
  indices[0] = u-1;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k =- 1; k>=delta_r; --k){
      index = arma::as_scalar(arma::find(rho == rho(u-1) + k));
      proposal(index) += 1;
      indices[-(k)] = index;
    }
  }


  return Rcpp::List::create(Rcpp::Named("proposal") = proposal,
                            Rcpp::Named("indices") = indices,
                            Rcpp::Named("delta_r") = delta_r,
                            Rcpp::Named("prob_forward") = prob_forward,
                            Rcpp::Named("prob_backward") = prob_backward
                              );
}

// Compute the distance between all rows in R and rho
double rank_dist_matrix(arma::mat R, arma::rowvec rho, std::string metric){
  int N = R.n_rows;

  if(R.n_cols != rho.n_elem) Rcpp::stop("R and rho have different number of elements");

  double total_distance = 0;

  for(int i = 0; i < N; ++i){
    total_distance += get_rank_distance(R.row(i), rho, metric = metric);
  }

  return total_distance;
}




