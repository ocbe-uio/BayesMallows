#include <RcppArmadillo.h>

double rtruncbeta(double shape1, double shape2, double trunc) {
  double Fa = R::pbeta(0, shape1, shape2, true, false);
  double Fb = R::pbeta(trunc, shape1, shape2, true, false);
  double u = R::runif(Fa, Fb);
  return R::qbeta(u, shape1, shape2, true, false);
}
