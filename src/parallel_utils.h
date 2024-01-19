#pragma once
#include <algorithm>

#ifdef __TBB_H
#include <tbb/parallel_for_each.h>
#endif

template <typename Iterator, typename Function>
void par_for_each(Iterator first, Iterator last, Function func) {
#ifdef __TBB_H
  tbb::parallel_for_each(first, last, func);
#else
  std::for_each(first, last, func);
#endif
}
