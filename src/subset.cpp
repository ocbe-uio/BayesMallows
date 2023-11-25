// SUBSET is a C++ library of Combinatorial Subroutines developed by John Burkardt, distributed under GPL.
// These functions are downloaded from http://people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html
// and have been modified to work with R and Rcpp.

# include <Rcpp.h>
# include <cstdlib>
# include <cstring>
# include <ctime>


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

int i4_sign ( int i )
{
  if ( i < 0 )
  {
    return (-1);
  }
  else
  {
    return 1;
  }
}

bool perm0_check ( int n, int p[] )
{
  bool check;
  int location;
  int value;

  check = true;

  for ( value = 0; value < n; value++ )
  {
    check = false;

    for ( location = 0; location < n; location++ )
    {
      if ( p[location] == value )
      {
        check = true;
        break;
      }
    }

    if ( ! check )
    {
      Rcpp::Rcout << "\n";
      Rcpp::Rcout << "PERM0_CHECK - Warning!\n";
      Rcpp::Rcout << "  Permutation is missing value " << value << "\n";
      break;
    }
  }

  return check;
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

void perm0_mul ( int n, int p1[], int p2[], int p3[] )
{
  int i;

  if ( !perm0_check ( n, p1 ) )
  {
    //cerr << "\n";
    //cerr << "PERM0_MUL - Fatal error!\n";
    //cerr << "  PERM0_CHECK rejects permutation.\n";
    Rcpp::stop("error");
  }

  if ( !perm0_check ( n, p2 ) )
  {
    //cerr << "\n";
    //cerr << "PERM0_MUL - Fatal error!\n";
    //cerr << "  PERM0_CHECK rejects permutation.\n";
    Rcpp::stop("error");
  }

  for ( i = 0; i < n; i++ )
  {
    p3[i] = p2[p1[i]];
  }

  return;
}

int *perm0_inverse ( int n, int p1[] )
{
  int i;
  int i0;
  int i1;
  int i2;
  int is;
  int *p2;

  if ( n <= 0 )
  {
    Rcpp::stop("error");
  }

  if ( !perm0_check ( n, p1 ) )
  {
    Rcpp::stop("error");
  }

  p2 = new int[n];
  for ( i = 0; i < n; i++ )
  {
    p2[i] = p1[i] + 1;
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p2[i-1];

    while ( i < i1 )
    {
      i2 = p2[i1-1];
      p2[i1-1] = - i2;
      i1 = i2;
    }

    is = - i4_sign ( p2[i-1] );
    p2[i-1] = abs ( p2[i-1] ) * i4_sign ( is );
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

int perm0_distance ( int n, int a[], int b[] )
{
  int *binv;
  int *c;

  int length;
  int *sub;
  int value;

  c = new int[n];
  sub = new int[n];

  binv = perm0_inverse ( n, b );

  perm0_mul ( n, a, binv, c );

  perm_ascend ( n, c, length, sub );

  delete [] binv;
  delete [] c;
  delete [] sub;

  value = n - length;

  return value;
}
