// SUBSET is a C++ library of Combinatorial Subroutines developed by John Burkardt, distributed under GPL.
// These functions are downloaded from http://people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html
// and have been modified considerably to use RcppArmadillo rather than C style arrays.

# include <RcppArmadillo.h>
using namespace arma;

int perm_ascend ( const ivec& a) {
  int length{};
  int n = a.size();
  ivec top(n, fill::zeros);
  ivec top_prev(n, fill::zeros);

  for (int i{}; i < n; i++) {
    int k = 0;
    for (int j{}; j < length; j++) {
      if (a(i) <= a(top(j) - 1)) {
        k = j + 1;
        break;
      }
    }

    if (k == 0) {
      length += 1;
      k = length;
    }
    top(k - 1) = i + 1;
    top_prev(i) = k > 1 ? top(k - 2) : 0;
  }

  return length;
}

ivec perm0_mul ( const ivec& p1, const ivec& p2) {
  const unsigned int n = p1.size();
  ivec p3(n);
  for (size_t i{}; i < n; i++) p3(i) = p2(p1(i));
  return p3;
}

ivec perm0_inverse (const ivec& p1) {
  int n = p1.size();
  ivec p2 = p1 + 1;

  for (int i{}; i < n; i++) {
    int i1 = p2(i);
    while ( (i + 1) < i1 ) {
      int i2 = p2(i1 - 1);
      p2(i1-1) = - i2;
      i1 = i2;
    }
    p2(i) = abs(p2(i)) * ((- p2(i) < 0) ? -1 : 1);
  }

  for (int i{}; i < n; i++) {
    int i1 = - p2(i);
    if (0 <= i1) {
      int i0 = i + 1;
      for ( ; ; ) {
        int i2 = p2(i1 - 1);
        p2(i1 - 1) = i0;
        if (i2 < 0) break;
        i0 = i1;
        i1 = i2;
      }
    }
  }
  return p2 - 1;
}

int perm0_distance ( const ivec& a, const ivec& b ) {
  int n = a.size();
  ivec binv = perm0_inverse ( b );
  ivec c = perm0_mul (a, binv);
  int length = perm_ascend(c);

  return n - length;
}
