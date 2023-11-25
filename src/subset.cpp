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

void i4vec_decrement ( int n, int v[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = v[i] - 1;
  }

  return;
}

void perm_ascend ( int n, int a[], int &length, int sub[] )
{
  int i;
  int j;
  int k;
  int *top;
  int *top_prev;

  if ( n <= 0 )
  {
    length = 0;
    return;
  }

  top = new int[n];
  for ( i = 0; i < n; i++ )
  {
    top[i] = 0;
  }

  top_prev = new int[n];
  for ( i = 0; i < n; i++ )
  {
    top_prev[i] = 0;
  }
  for ( i = 0; i < n; i++ )
  {
    sub[i] = 0;
  }

  length = 0;

  for ( i = 1; i <= n; i++ )
  {
    k = 0;

    for ( j = 1; j <= length; j++ )
    {
      if ( a[i-1] <= a[top[j-1]-1] )
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
      top_prev[i-1] = top[k-2];
    }
    else
    {
      top_prev[i-1] = 0;
    }
  }

  j = top[length-1];
  sub[length-1] = a[j-1];

  for ( i = length-1; 1 <= i; i-- )
  {
    j = top_prev[j-1];
    sub[i-1] = a[j-1];
  }

  delete [] top;
  delete [] top_prev;

  return;
}

void perm0_mul ( const arma::ivec& p1, int p2[], int p3[] )
{
  int i;
  int n = p1.size();
  for ( i = 0; i < n; i++ )
  {
    p3[i] = p2[p1(i)];
  }

  return;
}

int *perm0_inverse ( const arma::ivec& p1 )
{
  int i;
  int i0;
  int i1;
  int i2;
  int *p2;
  int n = p1.size();

  p2 = new int[n];
  for ( i = 0; i < n; i++ )
  {
    p2[i] = p1(i) + 1;
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = p2[i-1];

    while ( i < i1 )
    {
      i2 = p2[i1-1];
      p2[i1-1] = - i2;
      i1 = i2;
    }

    p2[i-1] = abs ( p2[i-1] ) * ((- p2[i-1] < 0) ? -1 : 1);
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p2[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p2[i1-1];
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

  i4vec_decrement ( n, p2 );

  return p2;
}

int perm0_distance ( const arma::ivec& a, const arma::ivec& b )
{
  int *binv;
  int *c;

  int length;
  int *sub;
  int value;

  int n = a.size();
  c = new int[n];
  sub = new int[n];

  binv = perm0_inverse ( b );

  perm0_mul ( a, binv, c );

  perm_ascend ( n, c, length, sub );

  delete [] binv;
  delete [] c;
  delete [] sub;

  value = n - length;

  return value;
}
