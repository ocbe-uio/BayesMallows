#include "classes.h"
#include "distances.h"
using namespace arma;

Clustering::Clustering(const Parameters& pars,
                       const Rcpp::List& compute_options,
                       const unsigned int n_assessors) :
  clustering {pars.n_clusters > 1},
  clus_thinning { compute_options["clus_thinning"] },
  index { regspace<uvec>(0, pars.n_clusters - 1) },
  include_wcd { compute_options["include_wcd"] },
  save_ind_clus { compute_options["save_ind_clus"] } {
    int n_cluster_assignments = pars.n_clusters > 1 ?
    std::ceil(pars.nmc * 1.0 / clus_thinning) : 1;
    cluster_probs = zeros(pars.n_clusters, n_cluster_assignments);
    cluster_probs.col(0).fill(1.0 / pars.n_clusters);
    current_cluster_probs = cluster_probs.col(0);
    cluster_assignment = zeros<umat>(n_assessors, n_cluster_assignments);
    ivec a = Rcpp::sample(pars.n_clusters, n_assessors, true) - 1;
    cluster_assignment.col(0) = conv_to<uvec>::from(a);
    current_cluster_assignment = cluster_assignment.col(0);

    dist_mat = zeros(n_assessors, pars.n_clusters);
    within_cluster_distance = zeros(
      pars.n_clusters, include_wcd ? pars.nmc : 1);
    update_wcd(0);
  }

void Clustering::update_cluster_probs(
    const Parameters& pars,
    const Priors& pris){
  vec cluster_probs(pars.n_clusters);

  for(size_t i{}; i < pars.n_clusters; ++i){
    cluster_probs(i) = R::rgamma(
      sum(current_cluster_assignment == i) + pris.psi, 1.0);
  }
  current_cluster_probs = normalise(cluster_probs, 1);
}

void Clustering::update_cluster_labels(
    const int t,
    const Data& dat,
    const Parameters& pars,
    const std::unique_ptr<PartitionFunction>& pfun
){
  uvec new_cluster_assignment(dat.n_assessors);
  mat assignment_prob(dat.n_assessors, pars.n_clusters);
  for(size_t i{}; i < pars.n_clusters; ++i){
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      pars.alpha_old(i) / dat.n_items * dist_mat.col(i) -
      pfun->logz(pars.alpha_old(i));
  }

  for(size_t i = 0; i < dat.n_assessors; ++i){
    rowvec probs = exp(assignment_prob.row(i) - max(assignment_prob.row(i)) -
      log(sum(exp(assignment_prob.row(i) - max(assignment_prob.row(i))))));

    ivec ans(probs.size());
    R::rmultinom(1, probs.begin(), probs.size(), ans.begin());
    new_cluster_assignment(span(i)) = find(ans == 1);
  }

  if(save_ind_clus){
    assignment_prob.save(
      std::string("cluster_probs") + std::to_string(t + 1) +
        std::string(".csv"), csv_ascii);
  }
  current_cluster_assignment = new_cluster_assignment;
}

void Clustering::update_wcd(const int t){
  if(!include_wcd) return;

  const unsigned int n_clusters = dist_mat.n_cols;
  vec wcd(n_clusters);

  for(size_t i{}; i < n_clusters; ++i){
    mat dist_vec = dist_mat.submat(
      find(current_cluster_assignment == i), index.subvec(i, i));
    wcd(i) = accu(dist_vec);
  }
  within_cluster_distance.col(t) = wcd;
}

void Clustering::update_dist_mat(
    const Data& dat, const Parameters& pars,
    const std::unique_ptr<Distance>& distfun){
  if(clustering | include_wcd) {
    for(size_t i{}; i < pars.n_clusters; ++i)
      dist_mat.col(i) = distfun->matdist(dat.rankings, pars.rho_old.col(i));
  }
}
