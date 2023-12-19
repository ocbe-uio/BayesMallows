// #pragma once
// #include <memory>
// #include <RcppArmadillo.h>
// #include <Rmath.h>
// #include "distances.h"
//
// struct AugmentationProposal {
//   AugmentationProposal() {};
//   virtual ~AugmentationProposal() = default;
//   virtual void propose_augmentation() = 0;
//
//   double log_hastings_correction{};
//   arma::vec proposal{};
//
//   arma::uvec unranked_items{};
//
//   arma::uvec sample_inds(const arma::uvec& unranked_items) {
//     Rcpp::IntegerVector inds =
//       Rcpp::sample(unranked_items.size(), unranked_items.size()) - 1;
//     return Rcpp::as<arma::uvec>(Rcpp::wrap(inds));
//   }
// };
//
//
//
// struct UniformProposal : AugmentationProposal {
//   void propose_augmentation(
//       const arma::vec& ranks, const arma::uvec& indicator) {
//     unranked_items = arma::find(indicator == 1);
//     arma::uvec inds = sample_inds(unranked_items);
//     proposal = ranks;
//
//     arma::vec available_rankings = ranks(unranked_items);
//
//     proposal(unranked_items) = available_rankings(inds);
//   }
// };
//
// struct PseudoLikelihoodProposal : AugmentationProposal {
//   void propose_augmentation(
//       const arma::vec& ranks, const arma::uvec& indicator,
//       const double alpha, const arma::vec& rho,
//       const std::unique_ptr<Distance>& distfun) {
//     unranked_items = arma::find(indicator == 1);
//     arma::uvec inds = sample_inds(unranked_items);
//     unranked_items = unranked_items(inds);
//     proposal = ranks;
//
//     while(unranked_items.n_elem > 0) {
//       int item_to_rank = unranked_items(0);
//       arma::vec available_rankings = proposal(unranked_items);
//       arma::vec log_numerator(available_rankings.n_elem);
//       for(size_t ll{}; ll < available_rankings.n_elem; ll++) {
//         const arma::vec& v = available_rankings(arma::span(ll));
//         log_numerator(ll) = -alpha / ranks.n_elem *
//           distfun->d(v, rho(arma::span(item_to_rank)));
//       }
//       arma::vec sample_probs = arma::normalise(arma::exp(log_numerator), 1);
//
//       arma::ivec ans(sample_probs.size());
//       R::rmultinom(1, sample_probs.begin(), sample_probs.size(), ans.begin());
//       proposal(arma::span(item_to_rank)) =
//         available_rankings(arma::find(ans == 1));
//     }
//
//
//
//   }
// };
