// SUBSET is a C++ library of Combinatorial Subroutines developed by John Burkardt, distributed under GPL.
// These functions are downloaded from http://people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html
// and have been modified to work with R and Rcpp.

# include <RcppArmadillo.h>
# include <cstdlib>
# include <cstring>
# include <ctime>
#include "setdiff.h"


using namespace std;

# include "subset.h"

int perm_ascend ( const arma::ivec& a)
{
  int length{};
  int n = a.size();

  if ( n <= 0 ) return length;

  arma::ivec top(n, arma::fill::zeros);
  arma::ivec top_prev(n, arma::fill::zeros);

  for ( int i{1}; i <= n; i++ )
  {
    int k = 0;
    for ( int j{1}; j <= length; j++ )
    {
      if ( a(i-1) <= a(top(j-1)-1) )
      {
        k = j;
        break;
      }
    }

    if ( k == 0 )
    {
      length = length + 1;
      k = length;
    }

    top[k-1] = i;

    if ( 1 < k )
    {
      top_prev(i-1) = top(k-2);
    }
    else
    {
      top_prev(i-1) = 0;
    }
  }

  int j = top(length-1);

  for ( int i{length-1}; 1 <= i; i-- )
  {
    j = top_prev(j-1);
  }

  return length;
}

arma::ivec perm0_mul ( const arma::ivec& p1, const arma::ivec& p2)
{
  int n = p1.size();
  arma::ivec p3(n);
  for ( size_t i{}; i < n; i++ ) p3(i) = p2(p1(i));
  return p3;
}

arma::ivec perm0_inverse ( const arma::ivec& p1 )
{
  int i;
  int i0;
  int i1;
  int i2;
  int n = p1.size();

  arma::ivec p2 = p1 + 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p2(i-1);

    while ( i < i1 )
    {
      i2 = p2(i1-1);
      p2(i1-1) = - i2;
      i1 = i2;
    }

    p2(i-1) = abs ( p2(i-1) ) * ((- p2(i-1) < 0) ? -1 : 1);
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p2(i-1);

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p2(i1-1);
        p2[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }

        i0 = i1;
        i1 = i2;
      }
    }
  }

  return p2 - 1;
}

int perm0_distance ( const arma::ivec& a, const arma::ivec& b )
{

  int n = a.size();
  arma::ivec binv = perm0_inverse ( b );
  arma::ivec c = perm0_mul ( a, binv);
  int length = perm_ascend ( c);

  return n - length;
}
