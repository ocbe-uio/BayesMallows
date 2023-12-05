#include "classes.h"
using namespace arma;

Priors::Priors(
  const Rcpp::List& priors
) : lambda { verify_positive(Rcpp::as<double>(priors["lambda"])) },
kappa { Rcpp::as<ivec>(priors["kappa"]) },
psi { Rcpp::as<unsigned int>(priors["psi"]) } {}
