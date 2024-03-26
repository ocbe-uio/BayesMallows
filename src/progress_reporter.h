#pragma once

struct ProgressReporter{
  ProgressReporter(const Rcpp::List& progress_report);
  void report(size_t t);

private:
  const bool verbose;
  const size_t interval;
};
