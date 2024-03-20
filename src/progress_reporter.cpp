#include <Rcpp.h>
#include <cstddef>
#include "progress_reporter.h"

ProgressReporter::ProgressReporter(const Rcpp::List& progress_report) :
  verbose { progress_report["verbose"] },
  interval { progress_report["report_interval"] } {}

void ProgressReporter::report(size_t t) {
  if (t % interval == 0) {
    Rcpp::checkUserInterrupt();
    if(verbose){
      Rcpp::Rcout << "First " << t <<
        " iterations of Metropolis-Hastings algorithm completed." << std::endl;
    }
  }
}
