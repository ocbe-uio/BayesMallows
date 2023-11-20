#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "misc.h"

struct parameters {
  parameters(Rcpp::List model, Rcpp::List compute_options,
             Rcpp::List priors, Rcpp::List initial_values, int n_items) :
  n_items { n_items } {
    int nmc = compute_options["nmc"];
    int n_clusters = model["n_clusters"];
    int alpha_jump = compute_options["alpha_jump"];
    alpha.set_size(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
    double alpha_init = initial_values["alpha_init"];
    alpha.col(0).fill(alpha_init);
    alpha_old = alpha.col(0);

    int rho_thinning = compute_options["rho_thinning"];
    rho.set_size(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
    Rcpp::Nullable<arma::mat> rho_init = initial_values["rho_init"];
    rho.slice(0) = initialize_rho(n_items, n_clusters, rho_init);
    rho_old = rho(arma::span::all, arma::span::all, arma::span(0));

    kappa_1 = Rcpp::as<int>(priors["kappa_1"]);
    kappa_2 = Rcpp::as<int>(priors["kappa_2"]);

    error_model = Rcpp::as<std::string>(model["error_model"]);
    if(error_model == "bernoulli"){
      theta = arma::zeros<arma::vec>(nmc);
      shape_1 = arma::zeros<arma::vec>(nmc);
      shape_2 = arma::zeros<arma::vec>(nmc);
      shape_1(0) = kappa_1;
      shape_2(0) = kappa_2;
    } else {
      theta.reset();
      shape_1.reset();
      shape_2.reset();
    }


  }
  ~parameters() = default;

  void update_shape(int t, const arma::mat& rankings,
                    const Rcpp::List& constraints) {


    int n_assessors = rankings.n_cols;
    int sum_1{};
    int sum_2{};
    for(int i = 0; i < n_assessors; ++i){
      Rcpp::List assessor_constraints = Rcpp::as<Rcpp::List>(constraints[i]);
      for(int j = 0; j < n_items; ++j) {
        arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[j]);
        arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[j]);

        for(unsigned int k = 0; k < items_above.n_elem; ++k){
          int g = (arma::as_scalar(rankings.col(i).row(j)) < arma::as_scalar(rankings.col(i).row(items_above(k) - 1)));
          sum_1 += g;
          sum_2 += 1 - g;
        }
        for(unsigned int k = 0; k < items_below.n_elem; ++k){
          int g = (arma::as_scalar(rankings.col(i).row(j)) > arma::as_scalar(rankings.col(i).row(items_below(k) - 1)));
          sum_1 += g;
          sum_2 += 1 - g;
        }
      }
    }
    shape_1(t) = kappa_1 + sum_1;
    shape_2(t) = kappa_2 + sum_2;
    theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
  }

  arma::mat alpha;
  arma::vec alpha_old;
  arma::cube rho;
  arma::mat rho_old;
  int kappa_1;
  int kappa_2;
  std::string error_model;
  arma::vec theta;
  arma::vec shape_1;
  arma::vec shape_2;

private:
  int n_items;
};


#endif
