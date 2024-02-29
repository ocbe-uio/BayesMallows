#include <Rcpp.h>
#include <cstddef>
#include "progress_reporter.h"

ProgressReporter::ProgressReporter(bool verbose, size_t interval) :
  verbose { verbose }, interval { interval } {}

void ProgressReporter::report(size_t t) {
  if (t % 1000 == 0) {
    Rcpp::checkUserInterrupt();
    if(verbose){
      Rcpp::Rcout << "First " << t <<
        " iterations of Metropolis-Hastings algorithm completed." << std::endl;
    }
  }
}
