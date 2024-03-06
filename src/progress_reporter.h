#pragma once

struct ProgressReporter{
  ProgressReporter(bool verbose, size_t interval = 1000);
  void report(size_t t);

private:
  const bool verbose;
  const size_t interval;
};
